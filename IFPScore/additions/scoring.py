import json
import logging
import os.path
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from tqdm.auto import tqdm

from additions.descriptors import calc_plip_parallel
from additions.interactions import essential, interesting, required


def get_score(names, prefix="Descriptor_Custom_DataFrame_"):
    def calc(ifp):
        cols = [f"{prefix}{name}" for name in names]
        values = ifp[ifp.index.isin(cols)]
        return values.sum() / len(names)
    return calc


def get_ifp_scores(ifps, prefix="Descriptor_Custom_DataFrame_"):
    required_scores = 100 * ifps.apply(get_score(required, prefix=prefix), axis=1)
    essential_scores = 100 * ifps.apply(get_score(essential, prefix=prefix), axis=1)
    interesting_scores = 100 * ifps.apply(get_score(interesting, prefix=prefix), axis=1)
    return (1.5 * required_scores + 1.0 * essential_scores + 0.5 * interesting_scores) / 3.0, required_scores, essential_scores, interesting_scores


class DockingScorer:

    def __init__(self, docker, save_pose_info=False, dock_stereoisomers=True, output_dir=None):
        super().__init__()
        self.docker = docker
        self.save_pose_info = save_pose_info
        self.dock_stereoisomers = dock_stereoisomers
        if self.save_pose_info and not output_dir:
            raise ValueError("Output directory must be specified when saving pose info.")
        self.output_dir = output_dir
        self.counter = 0
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)

    def getScores(self, mols, frags=None):
        """
        Processes molecules and returns a score for each (i.e. a QSAR model prediction).
        """
        # gather info for docking and make some checks
        self.counter += 1
        pose_data = dict()
        info = dict()
        stereos_to_mols = dict()
        duplicates = dict()
        for idx, mol in tqdm(enumerate(mols), total=len(mols), desc="Preparing molecules for docking"):
            try:
                if type(mol) == str:
                    mol = Chem.MolFromSmiles(mol)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                smiles = Chem.MolToSmiles(mol)
                # prepare stereoisomers for docking
                stereoisomers = tuple(EnumerateStereoisomers(mol))
                stereoisomers_info = []
                if len(stereoisomers) > 1 and self.dock_stereoisomers:
                    for stereoisomer in stereoisomers:
                        stereoisomer = Chem.AddHs(stereoisomer)
                        AllChem.EmbedMolecule(stereoisomer)
                        stereoisomers_info.append((stereoisomer, Chem.MolToSmiles(stereoisomer, isomericSmiles=True)))
                if not stereoisomers_info:
                    stereoisomers_info = [(mol, smiles)]
                    if smiles not in stereos_to_mols:
                        stereos_to_mols[smiles] = {"smiles": smiles, "idx": idx}
                else:
                    for x in stereoisomers_info:
                        if x[1] not in stereos_to_mols:
                            stereos_to_mols[x[1]] = {"smiles": smiles, "idx": idx}
                # attach structures for docking
                for stereo_info in stereoisomers_info:
                    mol, smiles = stereo_info
                    if smiles not in info:
                        info[smiles] = {"idx": idx, "mol": mol}
                    else:
                        if smiles not in duplicates:
                            duplicates[smiles] = [idx]
                        else:
                            duplicates[smiles].append(idx)
                        raise ValueError(f"Duplicated smiles found: {smiles}, info entry exists: {info[smiles]}")
            except Exception as exp:
                logging.warning(exp, exc_info=True)
                info[str(idx)] = {"idx": idx}
                continue
        # dock molecules
        if len(info) > 0:
            mols_to_dock = [info[x]['mol'] for x in info if 'mol' in info[x]]
            results = self.docker.dock_molecule_list(mols_to_dock, 1)
            for result in results:
                if result[0] is not None:
                    pose = Chem.Mol(result[0], 0)
                    pose_smiles = Chem.MolToSmiles(pose, isomericSmiles=True)
                    if pose_smiles not in info:
                        logging.error(f"Pose yielded a different SMILES ({pose_smiles}) "
                                        f"from any of the inputs: {','.join(info.keys())}.")
                        continue
                    info[pose_smiles]["pose"] = pose
                    if self.save_pose_info:
                        pose_data[pose_smiles] = {
                            "metadata": result[1]
                        }
                        pose_data[pose_smiles]["original_mol"] = stereos_to_mols[pose_smiles]["smiles"]
                        pose_data[pose_smiles]["poses"] = []
                        for conf_id in range(result[0].GetNumConformers()):
                            pose_data[pose_smiles]["poses"].append(
                                Chem.MolToMolBlock(result[0], confId=conf_id)
                            )
        # calculate scores
        for conf_id in range(1):
            poses = [(smiles, Chem.Mol(values["pose"], confId=conf_id)) for smiles, values in info.items() if "pose" in values]
            if len(poses) != 0:
                ifps = calc_plip_parallel(poses, self.docker.protein.pdb_file, n_jobs=self.docker.n_cpus)
                scores, required_scores, essential_scores, interesting_scores = get_ifp_scores(ifps, prefix="")
                for smiles, score in scores.items():
                    if "score" not in info[smiles]:
                        info[smiles]["score"] = score
                    else:
                        info[smiles]["score"] = max(info[smiles]["score"], score)
                if self.save_pose_info:
                    ifps["SMILES"] = ifps.index.values
                    ifps["RequiredScores"] = required_scores
                    ifps["EssentialScores"] = essential_scores
                    ifps["InterestingScores"] = interesting_scores
                    ifps["TotalScores"] = scores
                    ifps["PosesAllData"] = [json.dumps(pose_data[smile]) for smile in ifps["SMILES"]]
                    ifps["Energy"] = [pose_data[smile]["metadata"][0]["vina_energy_total"] for smile in ifps["SMILES"]]
                    ifps.to_csv(os.path.join(self.output_dir, f"poses_{self.counter}.tsv"), sep="\t", index=False, header=True)
            else:
                logging.warning("Scoring ended without successful docking of any molecules.")
        # return scores
        scores = [0.0] * len(mols)
        for smiles in info.keys():
            if "score" in info[smiles]:
                idx = stereos_to_mols[smiles]["idx"]
                score = info[smiles]["score"]
                if scores[idx] < score:
                    scores[idx] = score
                    if smiles in duplicates:
                        for i in duplicates[smiles]:
                            scores[i] = score
        return np.array(scores)

    def getKey(self):
        """
        Unique Identifier among all the scoring functions used in a single environment.
        """

        return "DockingScorer"


class RemoteScriptScore:

    def __init__(self, name, nodes: list[str], n_cpus: list[int]):
        super().__init__()
        self.nodes = nodes
        self.n_cpus = n_cpus
        assert len(nodes) == len(n_cpus), "Number of nodes and number of CPUs must be the same."
        self.name = name

    def getScores(self, mols, frags=None):
        # divide mols to equal parts based on number of nodes
        mols = [Chem.MolToSmiles(mol, isomericSmiles=True) if type(mol) != str else mol for mol in mols]
        parts = np.array_split(mols, len(self.nodes))

        # send parts to nodes
        outfiles = []
        smi_files = []
        for part, node, n_cpus in zip(parts, self.nodes, self.n_cpus):
            name = f"{self.name}_{node}"

            # write smiles to file
            smiles_file = f"./data/temp/drugex/temp/docking_scores/{name}_smiles.txt"
            os.makedirs(os.path.dirname(smiles_file), exist_ok=True)
            if os.path.exists(smiles_file):
                os.remove(smiles_file)
            with open(smiles_file, "w") as f:
                f.write("\n".join(part))

            # run the docking script on the remote system
            scores_file = f"./data/temp/drugex/temp/docking_scores/{name}_scores.txt"
            if os.path.exists(scores_file):
                os.remove(scores_file)
            outfiles.append(scores_file)
            smi_files.append(smiles_file)
            command = f"cd ~/projects/drugex_docking/CCR5-drugex/ && screen -d -m nohup ./score_mols.sh {name} {n_cpus} &"
            os.system(f"ssh {node} -f '{command}'")

        # wait till all nodes finish
        for outfile in outfiles:
            while not os.path.exists(outfile):
                time.sleep(5)  # wait 5 seconds

        # collect results from file
        scores = []
        for smi_file, outfile in zip(smi_files, outfiles):
            with open(outfile, "r") as f:
                scores.extend([float(x) for x in f.readlines()])
            os.remove(outfile)
            os.remove(smi_file)

        return np.array(scores)

    def getKey(self):
        return "ScriptScore"
    
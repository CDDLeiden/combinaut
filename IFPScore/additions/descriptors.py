import os
import random
import tempfile
from typing import Iterable

import pandas as pd
from plip.structure.preparation import PDBComplex
from rdkit import Chem
from spock.parallel.iterators import parallel_iter, batched_iter
from tqdm.auto import tqdm


def read_pdb(pdb_file):
    # pdb_file = open(pdb_file, "r").read()
    # prot = Chem.MolFromMolBlock(pdb_file, sanitize=False, removeHs=False)
    prot = Chem.MolFromPDBFile(
        pdb_file,
        removeHs=False,
        sanitize=False,
        proximityBonding=True)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETAROMATICITY)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETCONJUGATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETHYBRIDIZATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SYMMRINGS)
    # Chem.SanitizeMol(prot, Chem.SANITIZE_PROPERTIES)
    Chem.SanitizeMol(prot, Chem.SANITIZE_CLEANUP)
    # prot = Chem.AddHs(prot, addCoords=True)
    return prot


def calc_plip(mols: Iterable[tuple[str, Chem.Mol]], pdb_file: str):
    prot = read_pdb(pdb_file)
    df = pd.DataFrame()
    for pose_id, mol in mols:
        complex = Chem.CombineMols(prot, mol)
        hash = random.getrandbits(128)
        temp_file = "%032x_complex.pdb" % hash
        temp_file = os.path.join(tempfile.gettempdir(), temp_file)
        Chem.MolToPDBFile(complex, temp_file, flavor=4)
        mol = PDBComplex()
        mol.load_pdb(temp_file)
        mol.analyze()
        longnames = [x.longname for x in mol.ligands]
        bsids = [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands]
        indices = [j for j, x in enumerate(longnames) if x == 'UNL']
        bsid = bsids[indices[0]]
        interactions = mol.interaction_sets[bsid].all_itypes
        df_mol = pd.DataFrame()
        for interaction in interactions:
            name = interaction.__class__.__name__.replace('_interaction', '')
            if hasattr(interaction,  'protisdon'):
                donor = "d" if interaction.protisdon else "a"
            else:
                donor = ''
            res_name = f'{name}{donor}_{interaction.restype}_{interaction.resnr}_{interaction.reschain}'
            df_mol[res_name] = True
        df_mol.loc[0, df_mol.columns] = True
        df_mol["Pose_ID"] = pose_id
        df = pd.concat([df, df_mol])
        os.remove(temp_file)

    df.fillna(False, inplace=True)
    return df


def calc_plip_parallel(mols: Iterable[tuple[str, Chem.Mol]], pdb_file: str, n_jobs: int = None, batch_size: int = 10):
    n_jobs = n_jobs if n_jobs else os.cpu_count()
    df_result = pd.DataFrame()
    for result in tqdm(parallel_iter(batched_iter(mols, batch_size), calc_plip, n_jobs, pdb_file=pdb_file), desc="Calculating PLIP interactions"):
        df_result = pd.concat([df_result, result])
    df_result.fillna(False, inplace=True)
    df_result.set_index("Pose_ID", inplace=True, drop=True, verify_integrity=True)
    return df_result

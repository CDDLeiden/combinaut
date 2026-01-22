import json
import logging
import os
import shutil

import pandas as pd
from rdkit import Chem
from spock.storage import TabularStorage
from tqdm.auto import tqdm

import settings

# load the storage with docked molecules
store = TabularStorage(settings.STORAGE_FOLDER)
target = store.get_target(settings.PROTEIN_NAME)

# define a smarts pattern here to analyze only a subset of molecules rather than the whole storage
smarts_pattern = ''
store = store.smarts_search(smarts_pattern) if smarts_pattern else store

# save poses and other data
poses_dir = f'{settings.DOCKING_FOLDER}/{smarts_pattern}'
if os.path.exists(poses_dir):
    shutil.rmtree(poses_dir)
os.makedirs(poses_dir)
with open(os.path.join(poses_dir, 'protein.pdb'), "w", encoding="utf-8") as f:
    f.write(target.get_pdb())
print(f"Matched {len(store)} molecules...")
df = pd.read_table(f"{settings.DATA_ROOT}/{settings.DATASET_NAME}.tsv")
df_store = []
df_poses = []
visited = set()
for mol in tqdm(store, total=len(store), desc=f"Saving docking data to {poses_dir}..."):
    # write all poses for each compound to an SDF (for easy visualization  with pymol)
    try:
        poses = store.get_poses(mol.id, target, raw=True)
    except ValueError:
        logging.warning(f"Failed to get poses for {mol.id}")
        continue
    sdfs_dir = os.path.join(poses_dir, 'sdfs')
    os.makedirs(sdfs_dir, exist_ok=True)
    poses_file = os.path.join(sdfs_dir, f'{mol.id}_poses.sdf')
    with open(poses_file, "w", encoding="utf-8") as f:
        f.write(poses)
    suppl = Chem.SDMolSupplier(poses_file)
    # select the first pose for each compound and save it to the poses data frame
    for pose in suppl:
        pose_id = pose.GetProp("ID")
        if pose_id.endswith("_0"):
            df_poses.append(pd.DataFrame(
                {
                    "Store_ID": [mol.id],  # ID of the compound in the storage
                    "Pose_ID": [pose_id],  # ID of the pose in the storage
                    "Energy": [float(pose.GetProp("vina_energy_total"))],
                    "SMILES": [mol.smiles],  # standardized SMILES
                    "SMILES_repr": [pose.GetProp("SMILES_isomeric")],  # docked SMILES
                    "SDF": [Chem.MolToMolBlock(pose)]  # sdf  of the best  docked pose
                }
            ))
    # collect data from the original dataset
    try:
        cids = [json.loads(mol.metadata['metadata'])['CID']]
        cids.extend([x['CID'] for x in json.loads(mol.metadata['~metadata_other'])])
        cids = [cid for cid in cids if cid not in visited]
        selection = df[df['CID'].isin(cids)].copy()
        selection['Store_ID'] = mol.id  # makes sure we can connect df_poses and df_store
        selection['SMILES'] = mol.smiles
        visited.update(cids)
    except KeyError:
        logging.warning(f"Failed to get metadata for {mol.id}. No CID column found. Skipping...")
        selection = pd.DataFrame({"Store_ID": [mol.id], "SMILES:": [mol.smiles]})
    df_store.append(selection)
df_store = pd.concat(df_store)
df_poses = pd.concat(df_poses)
data_file = os.path.join(poses_dir, 'data.tsv')
df_store.to_csv(data_file, sep='\t', index=False)
poses_file = os.path.join(poses_dir, 'poses.tsv')
df_poses.to_csv(poses_file, sep='\t', index=False)
print(f"Saved poses to {poses_file} and associated data to {data_file}")

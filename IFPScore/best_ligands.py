import os
import shutil

import pandas as pd
from rdkit import Chem

import settings
from additions.data import MolTable
from additions.descriptors import read_pdb

# load the data set with IFPs
mt = MolTable(
    name=f'{settings.PROTEIN_NAME}_ifps',
    store_dir=settings.DOCKING_FOLDER,
    index_cols=['Pose_ID'],
    smiles_col=f"SMILES_repr",
)
print(len(mt))

best_ids = []
for pattern in settings.ALLOSTERIC_PATTERNS:
    result = mt.searchWithSMARTS([pattern])
    result_df = result.getDF().sort_values(["pchembl_value_Median"], ascending=False).iloc[0:5, :]
    result_df = result_df[result_df.accessions.str.contains("P41597")]
    best_ids.extend(list(result_df['Pose_ID']))
    for idx, data in result_df.iterrows():
        print(data[mt.smilesCol], data['Pose_ID'], data['pchembl_value_Median'])
# best_mols_str = " ".join([f"./data/docking/{settings.PROTEIN_NAME}/sdfs/{mol}_poses.sdf" for mol in best_ids])
# pymol_command = f"pymol {best_mols_str} ./data/docking/{settings.PROTEIN_NAME}/protein.pdb"
# print(pymol_command)

# find interaction fingerprints for the IDs
# result_df = mt.getDF()
# result_df = result_df[result_df["pchembl_value_Median_P41597"] >= 6.5]
# best_ids = list(result_df['Pose_ID'])
mt_best = mt.searchOnColumn("Pose_ID", best_ids, name=f"{mt.name}_best_molecules")
mt_best.save()
ifps = mt_best.getDescriptors()
valid_cols = ifps.columns[ifps.apply(lambda x: not all(x == 0), axis=0)]
ifps = ifps[valid_cols]
print(ifps.sum().sort_values(ascending=False))

# save ifp data for best molecules
df_sdfs = pd.read_table(f"{settings.DOCKING_FOLDER}/poses.tsv")
df_sdfs = df_sdfs[df_sdfs["Pose_ID"].isin(best_ids)]
out_dir = f"{settings.DOCKING_FOLDER}/complexes_best"
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.makedirs(out_dir, exist_ok=True)
df_sdfs.to_csv(f"{out_dir}/best_poses.tsv", sep="\t", index=False)
ifps.to_csv(f"{out_dir}/best_ifps.tsv", sep="\t", index=True)

# save the complexes as PDBs in one file (can be used to run PLIP on the server)
prot = read_pdb(f'{settings.PROTEIN_FOLDER}/{settings.PROTEIN_NAME}.pdb')
for idx, row in df_sdfs.iterrows():
    mol = Chem.MolFromMolBlock(row.SDF)
    complex = Chem.CombineMols(mol, prot)
    Chem.MolToPDBFile(complex, f"{out_dir}/{row.Pose_ID}.pdb", flavor=4)

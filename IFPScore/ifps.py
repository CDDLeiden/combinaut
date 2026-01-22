import pandas as pd
from qsprpred.data.data import MoleculeTable
from qsprpred.data.utils.descriptorcalculator import CustomDescriptorsCalculator
from qsprpred.data.utils.descriptorsets import DataFrameDescriptorSet
from rdkit import Chem

import settings
from additions.descriptors import calc_plip_parallel

smarts_pattern = ''  # again set a pattern for molecules you want to include
docking_results_dir = f'{settings.DOCKING_FOLDER}/{smarts_pattern}'
df = pd.read_table(f'{docking_results_dir}/poses.tsv')
mt = MoleculeTable(
    df=df,
    store_dir=docking_results_dir,
    name=f'{settings.PROTEIN_NAME}_ifps',
    index_cols=['Pose_ID'],
    smiles_col=f"SMILES_repr",
    overwrite=True
)
df = mt.getDF()
df_fps = calc_plip_parallel(
    ((pose_id, Chem.MolFromMolBlock(sdf)) for pose_id,sdf in zip(df.index, df.SDF)),
    f'{docking_results_dir}/protein.pdb',
    settings.N_CPUS
)

# add ifps to molecule table
calc = CustomDescriptorsCalculator(
    desc_sets=[DataFrameDescriptorSet(
        df=df_fps
    )]
)
mt.addCustomDescriptors(calc, recalculate=True)
ifps = mt.getDescriptors()
print(f"{ifps.shape[1]} IFPs calculated for {ifps.shape[0]} compounds:")
print(ifps.columns)
mt.save()
print(f"The resulting molecule table ({mt.name}) was saved to {mt.storeDir}.")

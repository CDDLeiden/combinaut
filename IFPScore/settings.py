import os

from spock.docking.vina.cpu_local import VinaProtein

DATASET_NAME = 'CCRs_HUMAN_ALL'
STORAGE_NAME = f'{DATASET_NAME}_allosteric_enumerated_stereo'
DATA_ROOT = f'./data/{STORAGE_NAME}/data'
STORAGE_FOLDER = f'{DATA_ROOT}/storage/{STORAGE_NAME}'
# SMARTS patterns for intracellular allosteric ligands
ALLOSTERIC_PATTERNS = [
    "OC(=O)[R][R]O[R][R]NS(=O)(=O)[R]",  # i.e some specific sulfonamides
]
PROTEIN_NAME = "5T1A_clean_mutations_reversed_withHs"
configs = {
    "5T1A_clean_mutations_reversed_withHs": {
        "center": [5.1, 28.0, 187.6],
        "box_size": [16.2, 17.8, 17.4]
    }
}
VINA_CONFIG = configs[PROTEIN_NAME]
N_CPUS = os.cpu_count()  # number of cpus to use for docking
EXHAUSTIVENESS = 8  # Vina exhaustiveness parameter
SEED = 42  # random seed for random operations
PROTEIN_FOLDER = f'./data/proteins'
PROTEIN = VinaProtein(
    PROTEIN_NAME,
    f'{PROTEIN_FOLDER}/{PROTEIN_NAME}.pdb',
    f'{PROTEIN_FOLDER}/{PROTEIN_NAME}.pdbqt'
)
DOCKING_FOLDER = f'{DATA_ROOT}/docking/{PROTEIN_NAME}'

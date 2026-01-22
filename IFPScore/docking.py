from spock.storage import TabularStorage
from spock.docking.vina.cpu_local import VinaDockingCPULocal

import settings


store = TabularStorage(settings.STORAGE_FOLDER)
docking = VinaDockingCPULocal(
    protein=settings.PROTEIN,
    n_cpus=settings.N_CPUS,
    box_spec=settings.VINA_CONFIG,
    embed_mols=True,  # set to False if conformers are already generated
    exhaustiveness=settings.EXHAUSTIVENESS,
    seed=settings.SEED,
)
docking.dock_storage(store, chunk_size=1, overwrite=True, save=True)
store.summary()

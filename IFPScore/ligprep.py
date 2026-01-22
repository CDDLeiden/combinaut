import settings
from spock.storage import TabularStorage
from spock.prep import Dimorphite


store_path = settings.STORAGE_FOLDER
store = TabularStorage(store_path)
ligprep = Dimorphite()
ligprep.process_store(store)
store.save(force=True)
store.summary()

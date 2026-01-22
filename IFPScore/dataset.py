import os
import shutil

import pandas as pd
from spock.storage import TabularStorage
from spock.storage.base import InchiIdentifier
from spock.utils.standardizers.papyrus import PapyrusStandardizer

from additions.data import PapyrusExtended, moltable_to_store
import settings

# Papyrus data set settings
name = settings.DATASET_NAME  # name of the data set
quality = "low"  # get all available data
version = "latest"  # get the latest version of Papyrus

# extract accession keys from the first column of the data set
df_gene_list = pd.read_table("./data/pantherGeneList.txt", header=None)
acc_keys = df_gene_list[0].apply(lambda x: x.split("|")[-1].replace("UniProtKB=", "")).tolist()

# fetch data from Papyrus for the given targets (may take a while)
papyrus = PapyrusExtended(
    data_dir="./data",  # where to store the data
    stereo=True,  # get stereochemistry
    plus_only=False,  # make sure to use the whole Papyrus
    descriptors=None,  # do not download descriptors
)
dataset = papyrus.getData(
    acc_keys,
    quality,
    name=name,
    use_existing=True,  # use existing data set if it was already compiled before
    output_dir=settings.DATA_ROOT,
)

# add Morgan fingerprints to the data set and save
# desc_calculator = MoleculeDescriptorsCalculator(desc_sets=[
#     FingerprintSet(fingerprint_type="MorganFP", radius=3, nBits=2048)
# ])
# dataset.addDescriptors(desc_calculator, recalculate=False)
# dataset.save()

# return matching records to find allosteric only ligands
if settings.ALLOSTERIC_PATTERNS:
    result = dataset.searchWithSMARTS(settings.ALLOSTERIC_PATTERNS, name=f"{dataset.name}_allosteric")
    print(f"Patterns matched {len(result)} records.")
else:
    result = dataset
print(result.getSummary())
result.save()

# save the result to a standardized storage to use in docking
storage_path = settings.STORAGE_FOLDER
if os.path.exists(storage_path):
    shutil.rmtree(storage_path)
storage = TabularStorage(
    storage_path,
    standardizer=PapyrusStandardizer(),
    identifier=InchiIdentifier()
)
moltable_to_store(result, storage,  add_cols=[
    "CID", "accession", "all_doc_ids", "QSPRID"
], enum_stereo=True)  # enumerate stereochemistry (not always annotated in Papyrus)
storage.save()
storage.summary()

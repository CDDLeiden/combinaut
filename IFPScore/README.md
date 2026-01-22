# Docking and Analysis of CCR Ligands from Papyrus -- Derivation of IFP-based Scoring Function

This repo implements a pipeline for docking and analysis of putative CCR ligands from Papyrus. Each step is implemented as a separate script, and the pipeline can be run from start to finish by running the scripts in the order listed below. Before you start, you should have an environment with Python 3.10 or newer set up and the following packages installed within it:

```bash
# install our toolbox
pip install git+https://github.com/martin-sicho/papyrus-scaffold-visualizer.git@v0.6.0
git clone git@gitlab.services.universiteitleiden.nl:cdd/spock.git
cd spock && pip install -e . && cd ..
pip install git+https://github.com/CDDLeiden/QSPRpred.git@v2.1.0

# install vina
conda install -c conda-forge gcc
export LD_LIBRARY_PATH=/home/<USERNAME>/.conda/envs/<CONDA_ENV_NAME>/lib:$LD_LIBRARY_PATH
pip install meeko lxml dimorphite_dl
conda install -c conda-forge vina

# install plip
git clone git@github.com:pharmai/plip.git && cd plip && pip install --no-dependencies . && cd ..
```

After that you can run the pipeline by executing the scripts in order, you are encouraged to read the code and modify it to suit your needs. The scripts are as follows:

1. `python dataset.py` - downloads the dataset from Papyrus and performs a SMARTS search for CCR allosteric ligands, it also creates a storage to store prepared ligands and docking results (make sure you apply some basic filters here to remove compounds that are obviously not intracellular allosteric ligands)
2. `python ligprep.py` - prepares the ligands for docking in the storage
3. `python docking.py` - performs docking of the prepared ligands in the storage
4. `python docking_analysis.py` - analyzes the docking results and creates some result files in the `./data/docking` directory
5. `python ifps.py` - calculates the IFPs for the analyzed ligands
6. `python clustering.py` - performs clustering of the ligands based on their IFPs
9. `python best_ligands.py` - extracts the best ligands for the target of interest and analyzes their interactions -> this can be used to derive a good IF{-based score that will correlate better with activity than the dockings scores on their own
7. `python search.py` - performs a search for ligands similar in binding to the reference ligands by scoring them with IFPs, the `search_analysis.ipynb` notebook can be used to analyze the results in more detail and judge the quality of the score
8. `python pymol_clusters.py` - can be used to sort clusters based on average score and generates a command to open them in PyMOL (you can install pymol in a separate environment: `conda install -c conda-forge pymol-open-source`)

import json

import pandas as pd
from scaffviz.depiction.plot import Plot
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import minkowski

import settings
from additions.clustering import cluster_mols
from additions.data import MolTable
from additions.manifold import Manifold


def manhattan_distance(u, v):
    return minkowski(u, v, p=1)


# load the data set with IFPs
mt = MolTable(
    name=f'{settings.PROTEIN_NAME}_ifps',
    store_dir=settings.DOCKING_FOLDER,
    index_cols=['Pose_ID'],
    smiles_col=f"SMILES_repr",
)

# add some columns from our data table to show in the plot below
df_data = pd.read_table(f'{settings.DOCKING_FOLDER}/data.tsv')
activities_all = []
activities_labelled = []
docs = []
accessions = []
label_accessions = ["P41597"]
labelled = []
for idx, row in mt.getDF().iterrows():
    data = df_data[df_data["Store_ID"] == row["Store_ID"]]
    data_targets = data.loc[data["accession"].isin(label_accessions), :]
    median_activity_labelled = data_targets["pchembl_value_Median"].median()
    median_activity_all = data["pchembl_value_Median"].median()
    activities_labelled.append(median_activity_labelled)
    activities_all.append(median_activity_all)
    docs.append(json.dumps(sorted(list(data["all_doc_ids"]))))
    accessions_data = sorted(list(data["accession"]))
    accessions.append(json.dumps(accessions_data))
    if (set(label_accessions) & set(accessions_data)) and median_activity_labelled >= 6.5:
        labelled.append("_".join(label_accessions) + " Active")
    else:
        labelled.append("Other")
mt.addProperty("pchembl_value_Median", activities_all)
labelled_activity_prop = "pchembl_value_Median_" + "_".join(label_accessions)
mt.addProperty(labelled_activity_prop, activities_labelled)
mt.addProperty("all_doc_ids", docs)
mt.addProperty("accessions", accessions)
mt.addProperty("ActivityLabel", labelled)
mt.save()

# cluster if not already done
recluster = True
if recluster or not mt.hasProperty('Cluster'):
    X = mt.getDescriptors()
    Y = mt.getDF().Energy < -10.0
    metrics = ['jaccard', manhattan_distance]
    color_thresholds = [0.8, 20]  # for each metric a different threshold
    methods = ['complete']  # linkage methods
    linkage_map = cluster_mols(
        X,
        Y.to_numpy(),
        methods=methods,
        metrics=metrics,
        color_thresholds=color_thresholds,
        outfile='clusters.png'
    )
    fl = fcluster(linkage_map['jaccard']['complete'], 45, criterion='maxclust')
    mt.addProperty('Cluster', [f"Cluster {i}" for i in fl])
    mt.save()

# choose a manifold for 2D embedding
n_components = 2
from sklearn.manifold import TSNE
tsne = TSNE(perplexity=15, random_state=settings.SEED, metric="jaccard")
# from sklearn.decomposition import PCA
# pca = PCA(n_components=n_components, random_state=42)
# from umap import UMAP
# umap = UMAP(n_components=n_components, n_neighbors=30, min_dist=0.25, metric='correlation', random_state=42, n_jobs=1)

# plot the ifp space and show clusters
print("Plotting...")
plt = Plot(Manifold(tsne))
plt.plot(
    mt,
    color_by="ActivityLabel",
    recalculate=True,
    card_data=[labelled_activity_prop, "pchembl_value_Median", "Energy", "Cluster", "all_doc_ids", "accessions"],
    title_data='Pose_ID'
)
mt.save()  # save the embedding coordinates to prevent re-computation

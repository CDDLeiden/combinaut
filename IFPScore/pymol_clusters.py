import settings
from additions.data import MolTable

# generate a command for pymol to show molecules from the given clusters
clusters = ["Cluster 1",]

# load the data set with IFPs
mt = MolTable(
    name=f'{settings.PROTEIN_NAME}_ifps',
    store_dir=settings.DOCKING_FOLDER,
    index_cols=['Pose_ID'],
    smiles_col=f"SMILES_repr",
)

# show best clusters
df_best_clusters = mt.getDF()[["Cluster", "Score"]].groupby("Cluster").max().sort_values("Score", ascending=False)
print("Clusters by score:")
print(df_best_clusters)

# search stuff if needed
# mt = mt.searchOnColumn("Cluster", ["P41597"])
mt = mt.searchOnColumn("Cluster", clusters, exact=True)

df_best_clusters = mt.getDF()[["Cluster", "Score"]].groupby("Cluster").max().sort_values("Score")
# best_clusters = df_best_clusters[df_best_clusters.Energy < 0.0].index.tolist()
best_clusters = df_best_clusters.index.tolist()
best_mols = mt.getDF()[mt.getDF().Cluster.isin(best_clusters)]["Store_ID"].tolist()
best_mols_str = " ".join(sorted([f"{settings.PROTEIN_FOLDER}/sdfs/{mol}_poses.sdf" for mol in best_mols]))
pymol_command = f"pymol {best_mols_str} {settings.PROTEIN_FOLDER}/protein.pdb"
print("copy this to the terminal to run pymol:")
print(pymol_command)

"""
Search for a pattern in a molecule table.
"""
from scaffviz.depiction.plot import Plot
from sklearn.manifold import TSNE

import settings
from additions.data import MolTable
from additions.manifold import Manifold
from additions.scoring import get_ifp_scores

# calculate the scores
mt = MolTable(
    name=f'{settings.PROTEIN_NAME}_ifps',
    store_dir=settings.DOCKING_FOLDER,
    index_cols=['Pose_ID'],
    smiles_col=f"SMILES_repr",
)
ifps = mt.getDescriptors()

# calculate scores and save to the table
scores, required_scores, essential_scores, interesting_scores = get_ifp_scores(ifps)
mt.addProperty("RequiredScore", required_scores.to_list())
mt.addProperty("EssentialScore", essential_scores.to_list())
mt.addProperty("InterestingScore", interesting_scores.to_list())
mt.addProperty("Score", scores.to_list())
mt.save()

tsne = TSNE(perplexity=30, random_state=42, metric="jaccard", n_jobs=24)
plt = Plot(Manifold(tsne))
plt.plot(
    mt,
    color_by="Score",
    # color_by="Cluster",
    # color_by="pchembl_value_Median_P41597",
    # color_by="EssentialScore",
    # x="Score",
    # y="NumHeavyAtoms",
    recalculate=False,
    card_data=[
        "pchembl_value_Median_P41597",
        "pchembl_value_Median",
        "Energy",
        "Score",
        "RequiredScore",
        "EssentialScore",
        "InterestingScore",
        # "all_doc_ids",
        # "accessions",
        "Cluster",
        "ActivityLabel",
    ],
    title_data='Pose_ID'
)

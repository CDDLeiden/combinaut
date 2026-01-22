import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import dendrogram, linkage


def cluster_mols(X, y=None, methods=('euclidean',), metrics=('single'), color_thresholds=None, outfile=None):
    if y is not None:
        labels = set(y)
        colors = list(mcolors.BASE_COLORS.keys())[0:len(labels)]
        label_c = {str(label) : color  for label, color in zip(labels, colors)}
    else:
        label_c = None

    n_row = len(metrics)
    n_col = len(methods)
    fig, axes = plt.subplots(n_row, n_col, figsize=(15*n_col, 15*n_row))
    linkage_map = {}
    for r, metr in enumerate(metrics):
        linkage_map[metr] = {}
        for c, meth in enumerate(methods):
            ax = axes.flat[r*n_col + c]
            # The default metric is euclidean - make sure it's printed properly in the plot-title
            if metr is None:
                metr = 'euclidean'
            print("{0} : {1}, {2} : {3}".format(r,metr, c, meth))
            Z = linkage(X, method=meth, metric=metr)
            dendrogram(Z, ax=ax, labels=y if y is not None else None, color_threshold=color_thresholds[r] if color_thresholds is not None else None) # ax=axes[i,j],
            ax.set_title("{0} linkage, {1} metric".format(meth, metr))
            ax.set_ylabel('Distance')
            ax.set_xlabel('Label')
            linkage_map[metr][meth] = Z

    # Color all observations according to their labels
    if y is not None:
        for ax in axes.flat:
            xlbls = ax.get_xmajorticklabels()
            for lbl in xlbls:
                lbl.set_color(label_c[lbl.get_text()])

    if outfile:
        plt.savefig(outfile)

    return linkage_map

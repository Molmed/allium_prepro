import umap
import pandas as pd
import matplotlib.pyplot as plt
from .subtype_thesaurus import SubtypeThesaurus


class BatchUmap():
    def _colormap(items):
        # Get palette from seaborn
        colormap = plt.cm.get_cmap('tab20')

        # For each item, assign a color
        class_to_color = {
            cls: colormap(i) for i, cls in enumerate(items)}

        return class_to_color

    def __init__(self,
                 prefix,
                 counts_file,
                 batches_file,
                 phenotype_file,
                 output_dir,
                 do_transform=False):
        self._prefix = prefix
        self._counts_file = counts_file
        self._batches_file = batches_file
        self._phenotype_file = phenotype_file
        self._output_dir = output_dir
        self._FONT_SIZE = 15
        self._FIG_SIZE = (8, 8)
        self._do_transform = do_transform

    def run(self):
        pheno_df = pd.read_csv(self._phenotype_file, index_col=0, sep=';')
        batches_df = pd.read_csv(self._batches_file, index_col=0)

        # Get unique values of the batch column
        unique_batches = batches_df['batch'].unique()
        n_components = unique_batches.size

        # join the two dataframes
        joined_df = pheno_df.join(batches_df, how='inner')
        features = pd.DataFrame(
            data=joined_df[['subtype', 'batch']],
            columns=['subtype', 'batch']).reset_index(drop=True)
        counts_df = pd.read_csv(self._counts_file, index_col=0)

        if self._do_transform:
            counts_df = counts_df.T

        mapper = umap.UMAP(n_components=n_components)
        data = mapper.fit_transform(counts_df)
        cols = ['UMAP_' + str(c+1) for c in range(n_components)]
        datadf = pd.DataFrame(data, columns=cols)
        finaldf = pd.concat([datadf, features], axis=1)
        finaldf.index = counts_df.index

        # Get colormap for subtypes
        colormap = BatchUmap._colormap(SubtypeThesaurus().allium_subtypes())

        # step factor =2 so we compare 1-2, 3-4 etc.
        for comp in range(1, n_components + 1, 2):
            plt.figure(figsize=self._FIG_SIZE)
            plt.xlabel('UMAP_{}'.format(comp), fontsize=self._FONT_SIZE)
            plt.ylabel('UMAP_{}'.format(comp + 1 ), fontsize=self._FONT_SIZE)
            plt.title('UMAP representation labeled by cytogenetic subtype',
                      fontsize=self._FONT_SIZE)
            clusterings = list(colormap.keys())
            gencolor = list(colormap.values())
            for clustering, coloring in zip(clusterings, gencolor):
                indicesToKeep = finaldf['subtype'] == clustering
                plt.scatter(finaldf.loc[indicesToKeep, 'UMAP_{}'.format(comp)],
                            finaldf.loc[indicesToKeep, 'UMAP_{}'.format(
                                comp + 1)],
                            c=coloring,
                            s=50)

            plt.legend(clusterings)
            plt.grid()
            plt.savefig(f'{self._output_dir}/{self._prefix}_umap_subtypes.png')

        # Get colormap
        colormap = BatchUmap._colormap(unique_batches)

        # step factor =2 so we compare 1-2, 3-4 etc.
        for comp in range(1, n_components + 1, 2):
            plt.figure(figsize=self._FIG_SIZE)
            plt.xlabel('UMAP_{}'.format(comp), fontsize=self._FONT_SIZE)
            plt.ylabel('UMAP_{}'.format(comp + 1), fontsize=self._FONT_SIZE)
            plt.tick_params('x', labelsize=30)
            plt.tick_params('y', labelsize=30)

            for clustering, coloring in zip(colormap.keys(),
                                            colormap.values()):
                indicesToKeep = finaldf['batch'] == clustering
                plt.scatter(finaldf.loc[indicesToKeep, 'UMAP_{}'.format(comp)],
                            finaldf.loc[indicesToKeep, 'UMAP_{}'.format(
                                comp + 1)],
                            c=coloring,
                            s=100)

            plt.legend(colormap.keys(), fontsize=self._FONT_SIZE)
            plt.axis('off')
            plt.savefig(f'{self._output_dir}/{self._prefix}_umap_batches.png')

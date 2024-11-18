import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
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
                 output_dir,
                 do_transform=False,
                 batches_file=None):
        self._prefix = prefix
        self._counts_file = counts_file
        self._batches_file = batches_file
        self._output_dir = output_dir
        self._FONT_SIZE = 15
        self._FIG_SIZE = (8, 8)
        self._do_transform = do_transform

    def run(self):
        # Load your dataset
        data = pd.read_csv(self._counts_file, index_col=0)

        if self._do_transform:
            data = data.T

        if not self._batches_file:
            # Create a dummy batch column
            batch_labels = pd.DataFrame(index=data.index)
            batch_labels['batch'] = np.random.choice(['All data'], data.shape[0])
            batch_labels = batch_labels['batch']
        else:
            batch_labels = pd.read_csv(self._batches_file, index_col=0)['batch']

        colormap = BatchUmap._colormap(batch_labels.unique())

        # Sort the data and batch labels by sample name
        data = data.sort_index()
        batch_labels = batch_labels.sort_index()

        # Normalize and scale the data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)

        # Apply UMAP for dimensionality reduction
        reducer = umap.UMAP(n_neighbors=15,
                            min_dist=0.1,
                            n_components=2,
                            random_state=42)
        umap_embedding = reducer.fit_transform(scaled_data)

        # Create a DataFrame for visualization
        umap_df = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])
        umap_df['Batch'] = batch_labels.values

        # Plot UMAP
        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            x='UMAP1',
            y='UMAP2',
            hue='Batch',
            palette=colormap,
            data=umap_df,
            s=50
        )
        plt.title('UMAP of Gene Expression Data Showing Batch Effects')
        plt.legend(title='Batch', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')
        plt.tight_layout()

        plt.savefig(f'{self._output_dir}/{self._prefix}_umap_batches.png')

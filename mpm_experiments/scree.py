import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

datasets = ['diedrich', 'heinaniemi', 'jude', 'krali', 'tran']
# datasets = ['heinaniemi']
data_path = '/home/mariya/Data/for_allium/allium'


def scree(filename, output_file):
    # Load GEX data (Example: CSV format)
    # Replace 'gex_data.csv' with your actual file
    df = pd.read_csv(filename, index_col=0)

    # Standardize the data (important for PCA)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df)

    # Perform PCA
    pca = PCA(n_components=10)
    pca.fit(scaled_data)

    # Explained variance ratio
    explained_variance = pca.explained_variance_ratio_

    # Scree plot
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(explained_variance) + 1), explained_variance, marker='o', linestyle='-')
    plt.xlabel('Principal Component')
    plt.ylabel('Explained Variance Ratio')
    plt.title('Scree Plot')
    plt.grid(True)

    # Save plot
    plt.savefig(output_file)


for dataset in datasets:
    # 1. Make scree plot before and after batch correction
    # 2. Make scree plot before and after normalization

    raw_file = before_file = f'{data_path}/{dataset}.counts.raw.csv'
    batch_corrected_file = f'{data_path}/{dataset}.tmp.counts.batch_corrected.csv'
    normalized_file = f'{data_path}/{dataset}.tmp.counts.norm.nolog.csv'

    if os.path.exists(batch_corrected_file):
        scree(raw_file,
              f'{data_path}/{dataset}_scree_before_batch_correction.png')
        scree(batch_corrected_file,
              f'{data_path}/{dataset}_scree_after_batch_correction.png')

    before_file = batch_corrected_file
    if not os.path.exists(before_file):
        before_file = raw_file

    scree(before_file, f'{data_path}/{dataset}_scree_before_norm.png')
    scree(normalized_file, f'{data_path}/{dataset}_scree_after_norm.png')

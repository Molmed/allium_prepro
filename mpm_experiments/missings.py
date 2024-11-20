
import pandas as pd

datasets = ['diedrich', 'heinaniemi', 'jude', 'krali', 'lilljebjorn', 'tran']

for DATASET_PREFIX in datasets:
    print("Processing ", DATASET_PREFIX)
    # Load the data
    data = pd.read_csv(f"/home/mariya/Data/allium/{DATASET_PREFIX}.missing_genes.csv", index_col=0)
    sigs = pd.read_csv("/home/mariya/Development/allium/data/signatures/signature_genes_v3.csv", index_col=0)

    ids = set(sigs['Gene ID']) & set(data.index)

    print("Missing genes in signatures:")
    # Print all rows in sigs where the Gene ID is in ids
    cols = ['Gene ID', 'feature_importance_mean']

    # Drop all colums in sigs except cols
    sigs = sigs[cols]
    print(sigs[sigs['Gene ID'].isin(ids)])
    print('\n\n\n')

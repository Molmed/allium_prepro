import pandas as pd
import os


class GexConcatenator():
    def __init__(self,
                 prefix,
                 raw_data_dir,
                 output_dir,
                 sample_name_extractor_func,
                 filename_filter_func=None,
                 separator='\t'):
        self._raw_data_dir = raw_data_dir
        self._files = os.listdir(raw_data_dir)
        self._output_file_path = f'{output_dir}/{prefix}.counts.raw.csv'
        self._sample_name_extractor_func = sample_name_extractor_func
        self._filename_filter_func = filename_filter_func
        self._separator = separator

    def concatenate(self):
        dfs = {}

        for filename in self._files:
            # Skip files that don't pass the filter
            if self._filename_filter_func and not self._filename_filter_func(filename):
                continue

            path = os.path.join(self._raw_data_dir, filename)

            # Extract sample name from filename
            sample = self._sample_name_extractor_func(filename)

            # Convert space separated text file to df, use first col as index
            sample_df = pd.read_csv(path,
                                    sep=self._separator,
                                    header=None,
                                    index_col=0)

            # Remove index name
            sample_df.index.name = None

            # Call the other column "count"
            sample_df.columns = ['count']

            # Create a new column with the sample name
            dfs[sample] = sample_df['count']

        # Concatenate all dataframes, using key as column name
        data = pd.concat(dfs, axis=1)

        # Sort columns by name
        data = data.sort_index(axis=1)

        # Drop all rows whose index starts with "__"
        data = data[~data.index.str.startswith("__")]

        # Write to file
        data.to_csv(self._output_file_path)

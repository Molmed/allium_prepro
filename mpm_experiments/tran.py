from lib.gex_concatenator import GexConcatenator

data_path = '/home/mariya/Data/tran'
raw_data_dir = f'{data_path}/GSE181157_RAW/'


def sample_name_extractor(x):
    return x.split('_', 1)[1].split('.')[0]


gc = GexConcatenator('tran',
                     raw_data_dir,
                     data_path,
                     sample_name_extractor)
gc.concatenate()

import pandas as pd
from Bio import SeqIO

def parser (file:str, format:str = 'fasta'):
    '''
    :param file: filepath for the file to be parsed
    :param format: format of the file to be parsed, supports 30+ formats (from SeqIO), defaults to fasta
    :return: a pandas dataframe containing the parsed sequences for future analysis
    '''

    try:
        with open(file) as f:
            pass
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file}")

    data = {}
    output = pd.DataFrame()

    records = list(SeqIO.parse(file, format))
    if not records:
        raise ValueError(f"No valid records found in file: {file}")
    for record in records:
        for key, value in record.__dict__.items():
            if key not in data:
                data[key] = []
            if hasattr(value, '_data'):
                value = str(value)
            data[key].append(value)

    for key, values in data.items():
        if values and len(values) == len(records):
            output[key.strip('_')] = values

    if 'seq' in output.columns:
        output['length'] = output['seq'].str.len()
    return output

#parser('synthetic_mtDNA_dataset.fasta', 'fasta').to_csv('out.csv')

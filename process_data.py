"""
Converts csv data into binary mutation matrix with rows as sample cells and columns as mutations. 
Mutations are identified through gene, chromosome, and amino acid change. 
"""
import pandas as pd

def get_binary_matrix():
    # read in csv, extract needed cols
    df = pd.read_csv("raw_data_1.csv")
    df = df[['sample_ID', 'gene', 'chr', 'amino_acid_change']]

    # clean up amino_acid_change column and remove "p." prefix
    df['amino_acid_change'] = df['amino_acid_change'].apply(lambda s: s[2:] if s[:2] == "p." else s)

    # create mutation identification column
    df['mutation'] = df['gene'] + "_" + df['chr'] + "_" + df['amino_acid_change']

    # clean up sample_ID column: remove unnecessary prefixes and suffixes and convert to int
    df['sample_ID'] = df['sample_ID'].apply(lambda s: int(s[4:-4]))

    # creates binary matrix of presence / absence in NumPy
    m = pd.crosstab(df['sample_ID'], df['mutation']).to_numpy()

    return m
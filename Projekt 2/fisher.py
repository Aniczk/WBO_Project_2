from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
import pandas as pd
import os
import numpy as np
from scipy.special import binom
import scipy.stats as stats

def count_hits(file):
    # reading the CSV file
    df = pd.read_csv(file)

    hits = 0
    mishits = 0
    for col in df.columns:
        if col != "Sequence ID":
            hits += df[col].sum()


    all = df.shape[0]*(df.shape[1]-1)
    mishits = all-hits

    return hits, mishits

def hits_table(array):
    df = pd.DataFrame(array, columns=["Hits", "Mishits"])
    df.index = ["File_1", "File_2"]
    df.loc['Column_Total'] = df.sum(axis=0)
    df.loc[:,'Row_Total'] = df.sum(axis=1)
    return df


def p(a, df): 
    # a = all proteins in 2 files
    K = df.iloc[0,2] # number of proteins in the first file
    n = df.iloc[1,2] # number of proteins in the second file
    N = df.iloc[2,2] # number of proteins in two files
    return (binom(K, a) *binom(N-K, n-a))/binom(N, n)

def fisher_implementation(df):
    a = df.iloc[0,0]
    p_observed = p(a,df)

    p_list=[]
    K = df.iloc[0,2] # number of proteins in the first file
    for i in range(K+1):
        if p(i,df)<=p_observed:
            p_list.append(p(i,df))     # append these probabilites to p_list only if <= p_observed
            
    fisher=np.sum(p_list) # the sum of this list corresponds to the p-value 

    return fisher

def fisher_scipy(df):
    ar = df.loc[['File_1','File_2'], ['Hits', 'Mishits']].values
    oddsratio, pvalue = stats.fisher_exact(ar)  
    return pvalue


# infile_1 = 'do_3/scan_pfam_extend.csv'
# infile_2 = 'do_3/scan_pfam_input_z2.csv'


def main():
    warnings.simplefilter("ignore")

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file_1", help="CSV file.")
    parser.add_argument("-j", "--input_file_2", help="CSV file.")
    args = vars(parser.parse_args())

    infile_1 = args["input_file_1"]
    infile_2 = args["input_file_2"]

    location = os.getcwd()
    input_file_1 = os.path.join(location,infile_1)
    input_file_2 = os.path.join(location,infile_2)

    hits_1, mishits_1 = count_hits(input_file_1)
    hits_2, mishits_2 = count_hits(input_file_2)

    ar = np.array([[hits_1, mishits_1],[hits_2, mishits_2]])
    df_hits = hits_table(ar)
    fisher_val_implementation = fisher_implementation(df_hits)
    fisher_val_scipy = fisher_scipy(df_hits)
    print("Result from the implementation of the Fisher test: ")
    print(fisher_val_implementation, '\n')
    print("Result from the scipy: ")
    print(fisher_val_scipy)
    


if __name__ == "__main__":
    main()
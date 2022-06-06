import subprocess
import os
import shutil
from Bio import SeqIO, SearchIO
import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings

def create_folder(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def delete_folder(dir):
    location = os.getcwd()

    # path
    path = os.path.join(location, dir)

    # removing directory
    if os.path.exists(path):
        shutil.rmtree(path, ignore_errors=False)


def crerate_small_files(input_file):
    create_folder("tmp")
    record_ids = []
    with open(input_file) as handle:
        i = 1
        for record in SeqIO.parse(handle, "fasta"):
            with open('tmp/fragment{}.fasta'.format(i), 'w') as f:
                f.write('>{}\n'.format(record.id))
                f.write('{}\n'.format(record.seq))
                f.close()
                record_ids.append(record.id)
                i+=1
    return record_ids

def hmmscan():
    create_folder("results")
    create_folder("results/hmmscan")

    path_of_the_directory= 'tmp'
    i = 1
    for filename in os.listdir(path_of_the_directory):
        input_file = os.path.join(path_of_the_directory,filename)
        if os.path.isfile(input_file):
            subprocess.run("curl -L -H 'Expect:' -H 'Accept:text/xml' -F hmmdb=pfam -F seq='<{}' https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan > 'results/hmmscan/scan_pfam_hmmscan{}.xml'".format(input_file,i), shell=True)
            i+=1

def create_csv_with_results(record_ids):
    """Use hmmscan outputs and create csv file with columns:
    Sequence ID, Domain1, Domain2, ...
    At the intersection of the row and column put
    0 if the domain was not found in the protein
    1 otherwise
    """

    list_of_files = os.listdir('results/hmmscan')
    df = {}
    for i in range(1,len(list_of_files)+1):
        f = open('results/hmmscan/scan_pfam_hmmscan{}.xml'.format(i),"r")
        lines = f.readlines()
        for line in lines:
            if " acc=" in line:
                domain = list(line.split(" "))[6][4:]
                if domain not in df:
                    L = [0]*len(list_of_files)
                    L[i-1]+=1
                    df[domain] = L
                else:
                    df[domain][i-1]+=1

    df["Sequence ID"] = record_ids
    data = pd.DataFrame(df)
    data = data.set_index("Sequence ID")
    # save results
    data.to_csv('results/scan_pfam.csv')


def main():
    warnings.simplefilter("ignore")

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--fasta_input_file", default="inputs/input-z2.fasta", help="Input file with protein sequences.")
    args = vars(parser.parse_args())

    infile = args["fasta_input_file"]

    location = os.getcwd()
    input_file = os.path.join(location,infile)
    rec_ids = crerate_small_files(input_file)
    hmmscan()
    delete_folder("tmp")
    create_csv_with_results(rec_ids)


if __name__ == "__main__":
    main()
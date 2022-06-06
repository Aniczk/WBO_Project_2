from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import os
import shutil
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings

def search_blast(file):
    """
    input: file - fasta file with sequences
    """
    create_folder("tmp")
    location = os.getcwd()

    outfile = os.path.join(location,"tmp","tmp_results.xml")

    sequence_data = open(file).read() 
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence_data)
    with open(outfile, 'w') as save_file: 
        blast_results = result_handle.read() 
        save_file.write(blast_results)
        save_file.close()

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

def find_similar_sequences(E_VALUE_THRESH, IDENTITY_TRESH):
    location = os.getcwd()

    blast_results_xml = os.path.join(location, "tmp","tmp_results.xml")

    records = NCBIXML.parse(open(blast_results_xml))

    sequences = []
    for r in records:
        for align in r.alignments:
            for hit in align.hsps:
                genome = align.accession
                organism = align.hit_def
                identity = hit.identities / hit.align_length
                if hit.expect <= E_VALUE_THRESH and identity >= IDENTITY_TRESH:
                    seq = hit.sbjct
                    sequences.append(seq)
    return sequences

def save_results(outfile, sequences):
    create_folder("results")
    outfile = os.path.join("results",outfile)

    with open(outfile, 'w') as f:
        for ID, seq in enumerate(sequences):
            f.write('>{}\n'.format(ID+1))
            f.write('{}\n'.format(seq))
        f.close()

def main():
    warnings.simplefilter("ignore")

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--fasta_input_file", default="inputs/input-z2-small.fasta", help="Input file with protein sequences.")
    parser.add_argument("-o", "--fasta_output_file", default="extend.fasta", help="Output fasta file.")
    parser.add_argument("-e", "--e_value_tresh", default=1e-10, help="e_value_tresh")
    parser.add_argument("-v", "--identity_value_tresh", default=1e-10, help="identity_value_tresh")

    args = vars(parser.parse_args())

    infile = args["fasta_input_file"]

    location = os.getcwd()
    input_file = os.path.join(location,infile)
    
    outfile = args["fasta_output_file"]
    E_VALUE_THRESH = args["e_value_tresh"]
    IDENTITY_TRESH = args["identity_value_tresh"]

    search_blast(input_file)
    sequences = find_similar_sequences(E_VALUE_THRESH, IDENTITY_TRESH)
    save_results(outfile, sequences)
    delete_folder("tmp")

if __name__ == "__main__":
    main()
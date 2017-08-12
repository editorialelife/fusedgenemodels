#!/usr/bin/python
from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import fileinput
import copy

import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

import argparse
def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Split proteins into C-terminal and N-terminal thirds ')
    parser.add_argument("-s", "--sequence", dest="sequence", help="fasta sequence file")
    args = parser.parse_args()

    if not args.sequence:
        parser.print_help()
        parser.error('Please specify a sequence file with -s or --sequence')

    
    try:
        with open(args.sequence,"rU") as seqFileHandle:
            proteinSequences = SeqIO.to_dict(SeqIO.parse(seqFileHandle, "fasta"))
            seqFileHandle.close()
    except IOError:
        print("Could not read file:",args.sequence)
        sys.exit()

    nterms = ""
    cterms = ""
    for record in proteinSequences:
        protLen = len(proteinSequences[record])

        if protLen > 100:
            shortname = record.split("|")[1]
            nterms = nterms + ">" + shortname + "\n" + proteinSequences[record][:protLen/3].seq + "\n"
            cterms = cterms + ">" + shortname + "\n" + proteinSequences[record][-1*(protLen/3):].seq + "\n"

    
    with open("Nterms.fasta", "w") as output_handle:
        print(nterms, file=output_handle)
        output_handle.close()

    with open("Cterms.fasta", "w") as output_handle:
        print(cterms, file=output_handle)
        output_handle.close()

if __name__ == "__main__":
    main(sys.argv[1:])
exit()

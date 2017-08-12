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

def most_common(lst):
    return max(set(lst), key=lst.count)

def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Compare blast hits for N-terminal and C-terminal portions of a protein')
    parser.add_argument("-C", "--cterm", dest="cterm", help="C-term blast output outfmt = 2")
    parser.add_argument("-N", "--nterm", dest="nterm", help="N-term blast output outfmt = 2")
    parser.add_argument("-O", "--orthomcl", dest="orthomcl", help="orthomcl output")
    args = parser.parse_args()

    if not args.cterm:
        parser.print_help()
        parser.error('Please speicify an output file for the C-terminal sequences')

    if not args.nterm:
        parser.print_help()
        parser.error('Please speicify an output file for the N-terminal sequences')

    #Parse blast output for Cterm
    Queries = {}
    try:
        with open(args.cterm,"rU") as seqFileHandle:
            for line in seqFileHandle:
                line = line.strip()
                if line[:6] == "Query=":
                    query = line[7:].strip()
                    Queries[query] = {'Cterm':[],'Nterm':[]}
                if "|" in line:
                    subject = line.split("|")[1].split(" ")[0]
                    Queries[query]['Cterm'].append(subject)
            seqFileHandle.close()
    except IOError:
        print("Could not read file:",args.cterm)
        sys.exit()

    #Parse blast output for Nterm
    try:
        with open(args.nterm,"rU") as seqFileHandle:
            for line in seqFileHandle:
                line = line.strip()
                if line[:6] == "Query=":
                    query = line[7:].strip()
                if "|" in line:
                    subject = line.split("|")[1].split(" ")[0]
                    Queries[query]['Nterm'].append(subject)
            seqFileHandle.close()
    except IOError:
        print("Could not read file:",args.nterm)
        sys.exit()

    #Build table of orthogroups
    orthogroup_dict = {}
    try:
        with open(args.orthomcl,"rU") as seqFileHandle:
            for line in seqFileHandle:
                line = line.strip().split(" ")
                orthogroup = line.pop(0)
                for entry in line:
                    orthogroup_dict[entry.split("|")[1]] = orthogroup
            seqFileHandle.close()
    except IOError:
        print("Could not read file:",args.orthomcl)
        sys.exit()

    


    outputstring = 'Query\tNterm_hits\tCterm_hits\tMin_hits\tOverlap\tN_orthogroup\tC_orthogroup\tSame\tFusionScore'
    for query in Queries:
        Cterm_set = set(Queries[query]['Cterm'])
        Nterm_set = set(Queries[query]['Nterm'])
        union = list(Nterm_set | Cterm_set)
        intersection = list(Nterm_set & Cterm_set)
        min_hits = min(len(Cterm_set),len(Nterm_set))
        if len(union) > 0:
            overlap = (len(intersection) * 100) / len(union)
        else:
            overlap = 100

        C_orthohits = []
        N_orthohits = []
            
        for Cterm_hit in Queries[query]['Cterm']:
            if Cterm_hit in orthogroup_dict:
                C_orthohits.append(orthogroup_dict[Cterm_hit])
        for Nterm_hit in Queries[query]['Nterm']:
            if Nterm_hit in orthogroup_dict:
                N_orthohits.append(orthogroup_dict[Nterm_hit])
            
        if len(C_orthohits) > 0:
            #C_orthogroup = C_orthohits
            C_orthogroup = most_common(C_orthohits[:5])
        else:
            C_orthogroup = "none"

        if len(N_orthohits) > 0:
            N_orthogroup = N_orthohits
            N_orthogroup = most_common(N_orthohits[:5])
        else:
            N_orthogroup = "none"

        fusion_score = 0
        orthogroups_agree = (N_orthogroup == C_orthogroup)
        if not orthogroups_agree:
            fusion_score += 5
            if N_orthogroup == "none" or C_orthogroup == "none":
                fusion_score -= 1

        if overlap < 1:
            if min_hits > 10:
                fusion_score += 5
            elif min_hits > 2:
                fusion_score += 4
            elif min_hits <= 2:
                fusion_score += 3

        elif overlap < 10:
            if min_hits > 10:
                fusion_score += 3
            elif min_hits > 2:
                fusion_score += 2
            elif min_hits <= 2:
                fusion_score += 1

        elif overlap > 50:
            fusion_score -= 3




            
        outputstring += '\n' + query
        outputstring += '\t' + str(Queries[query]['Nterm']) + '\t' + str(Queries[query]['Cterm'])
        outputstring += '\t' + str(min_hits) + '\t' + str(overlap)
        outputstring += '\t' + str(N_orthogroup) + '\t' + str(C_orthogroup) + '\t' + str(orthogroups_agree)
        outputstring += '\t' + str(fusion_score)
        
        
    print(outputstring)

if __name__ == "__main__":
    main(sys.argv[1:])
exit()

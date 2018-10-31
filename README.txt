

### This code is associated with the paper from Coradetti et al., "Functional genomics of lipid metabolism in the oleaginous yeast Rhodosporidium toruloides". eLife, 2018. http://dx.doi.org/10.7554/eLife.32110


Fused Gene Model Search Scripts
This is a collection of scripts and commands used to find gene models likely to be erroneous fusions of two neighboring genes in the genome of the oleaginous yeast Rhodosporidium toruloides (also known as Rhodotorula toruloides).  These scripts are provided as supplementary material for the scientific article "Functional genomics of lipid metabolism in the oleaginous yeast Rhodosporidium toruloides" submitted to the journal eLife in September 2017. They are offered here to clarify our approach in revising this particular genome (publicly available at http://genome.jgi.doe.gov/Rhoto_IFO0880_4/) and for any re-use or revision for which they may have utility.

Prerequisites 
1. A list of groups of orthologous genes output from OrthoMCL software with genes from the target genome and several reference genomes (orthogroups.txt) 
2. Fasta format file of protein sequences from the target genome with identifiers matching those in orthogroups.txt (RTO3.fasta)
3. Fasta format file of protein sequences from the reference genomes with identifiers matching those in orthogroups.txt (RTO3.fasta)
4. A working installation of NCBI BLAST+


Steps
1. Build BLAST database from reference protein sequences
2. Extract C-terminal and N-terminal 30% of amino acid sequences (SplitProteins.py)
3. Search C-terminal/N-terminal segments agains the reference database (blastp)
4. Compare sets of BLAST hits for overlap of orthogroups between matches for C-terminal and N-terminal sequences (CompareBlast.py) 


Example Command Sequence
#1. Build BLAST database
makeblastdb -in OrthoMCL_proteins.fasta -dbtype prot -out RefProteins -parse_seqids
#2. Extract C/N-terminal sequences
python SplitProteins.py -s RTO3.fasta
#3. Search C/N-terms against reference database
blastp -db RefProteins -evalue 1e-10 -num_alignments 0 -seg yes -outfmt 2 -query Cterms.fasta -out Cterm_hits_local &
blastp -db RefProteins -evalue 1e-10 -num_alignments 0 -seg yes -outfmt 2 -query Nterms.fasta -out Nterm_hits_local_test &
#4. Compare sets of BLAST hits
python CompareBlast.py -C Cterm_hits_local -N Nterm_hits_local -O orthogroups.txt > orthogroup_comparisons.txt


Output
orthogroup_comparisons.txt is a tab-delimited file with columns for the input gene identifier, gene identifiers in the reference genomes matched by the N and C termini, the number of genes matched by the N and C termini, the dominant orthologous group matched by the N and C termini, a TRUE/FALSE flag if those two dominant ortho-groups are the same and a 'Fusion Score' rating the likelihood that the input gene model was an erroneous fusion of two or more neighboring genes.  Genes with a score of 10 are very likely fused models, 0 or less very unlikely to be fused gene models.


Copyright (C) 2018 The Regents of the University of California All
rights reserved.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.


Disclaimer

NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY NOR THE REGENTS OF THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

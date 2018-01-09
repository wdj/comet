
NOTE: see https://code.ornl.gov/wjd/genomics_gpu/tree/master/tools
for the reference copy of the files in this directory

Preprocessing example
---------------------

export PATH="/lustre/atlas1/bif102/proj-shared/comet:$PATH"

mkdir 28M_permuted

shuf 28M_original/28M.txt > 28M_permuted/28M.txt

preprocess 28M_permuted/28M.txt 28M_permuted/28M.bin

allele_labels.sh 28M_permuted/28M.txt 28M_permuted/28M_allele_labels.txt

labels.sh 28M_permuted/28M.txt 28M_permuted/28M_labels.txt

line_indices 28M_permuted/28M.txt 28M_permuted/28M_line_indices.txt


CoMeT execution example
-----------------------








Postprocessing example
----------------------

postprocess_all.sh 2 28M_permuted/28M_allele_labels.txt \
  28M_permuted/28M_labels.txt \
  $(ls outs_full2_6000_2way_0_1_200/out_????.bin)

OPTIONAL:

ccc_validate_all.sh 2 28M_permuted/28M.txt \
  28M_permuted/28M_line_indices.txt \
  $(ls outs_full2_6000_2way_0_1_200/out_????.txt)

validate_all.sh 2 $(ls outs_full2_6000_2way_0_1_200/out_????.txt)


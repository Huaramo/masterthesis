#!/bin/bash
#This bash script calls the functions in the whole project.

#Call the programs computing expression correlation v.s. sequence identity.
#python ./plants/plants.py 2> ./plants/err > ./plants/py &
#python ./worms/worms.py 2> ./worms/err > ./worms/py &
#python ./mammals/mammals.py --fasta './gos/all_protein_seq_mr.fasta' --rnaseq ./mammals/mouse_rat_soellner.tsv --homologs './gos/ortho_inparalogs_mr.txt' './gos/ortho_orthologs_mr.txt' './gos/ortho_co_orthologs_mr.txt' --out ./mammals/summary_soellner.tsv 2> ./mammals/err1 > ./mammals/py1 &
#python ./mammals/mammals.py --fasta './gos/all_protein_seq_mr.fasta' --rnaseq ./mammals/mouse_rat_fushan.tsv --homologs './gos/ortho_inparalogs_mr.txt' './gos/ortho_orthologs_mr.txt' './gos/ortho_co_orthologs_mr.txt' --out ./mammals/summary_fushan.tsv 2> ./mammals/err2 > ./mammals/py2 &
#python ./mammals/mammals.py --fasta './gos/all_protein_seq_hm.fasta' --rnaseq ./mammals/human_mouse.tsv --homologs './gos/ortho_inparalogs_hm.txt' './gos/ortho_orthologs_hm.txt' './gos/ortho_co_orthologs_hm.txt' --out ./mammals/summary_hm.tsv 2> ./mammals/err3 > ./mammals/py3 &    

#------------------------------------------------------
#Preprocessing: Compute the tables of conditional probability and the tables for the computation of functional similarities:
python ./gos/prob_go.py 2> ./gos/err > ./gos/py &

#------------------------------------------------------
#Expression correlation v.s. semantic similarity.
#python ./gos/expr_go.py --rnaseq ./mammals/human_mouse.tsv --homologs './gos/ortho_inparalogs_hm.txt' './gos/ortho_orthologs_hm.txt' './gos/ortho_co_orthologs_hm.txt' --maps './gos/human_uniprot.tsv' './gos/mouse_uniprot.tsv' 2> ./gos/err1 > ./gos/py1 &
#python ./gos/expr_go.py --rnaseq ./mammals/mouse_rat_soellner.tsv --homologs './gos/ortho_inparalogs_mr.txt' './gos/ortho_orthologs_mr.txt' './gos/ortho_co_orthologs_mr.txt' --maps './gos/mouse_uniprot.tsv' './gos/rat_uniprot.tsv' 2> ./gos/err2 > ./gos/py2 &
#python ./gos/expr_go.py --rnaseq ./mammals/mouse_rat_fushan.tsv --homologs './gos/ortho_inparalogs_mr.txt' './gos/ortho_orthologs_mr.txt' './gos/ortho_co_orthologs_mr.txt' --maps './gos/mouse_uniprot.tsv' './gos/rat_uniprot.tsv' 2> ./gos/err3 > ./gos/py3 &


#------------------------------------------------------
#GO-term analyses:
#python ./gos/semsims.py --fasta './gos/all_protein_seq_hm.fasta' --homologs './gos/ortho_inparalogs_hm.txt' './gos/ortho_orthologs_hm.txt' './gos/ortho_co_orthologs_hm.txt' --maps './gos/human_uniprot.tsv' './gos/mouse_uniprot.tsv'  2> ./gos/err4 > ./gos/py4 &

#python ./gos/semsims.py --fasta './gos/all_protein_seq_mr.fasta' --homologs './gos/ortho_inparalogs_mr.txt' './gos/ortho_orthologs_mr.txt' './gos/ortho_co_orthologs_mr.txt' --maps './gos/mouse_uniprot.tsv' './gos/rat_uniprot.tsv' 2> ./gos/err5 > ./gos/py5 & 

#python ./gos/semsims.py --fasta './gos/all_protein_seq_cp.fasta' --homologs './gos/ortho_inparalogs_cp.txt' './gos/ortho_orthologs_cp.txt' './gos/ortho_co_orthologs_cp.txt'  --maps './gos/cerev_uniprot.tsv' './gos/pombe_uniprot.tsv' 2> ./gos/err6 > ./gos/py6 &















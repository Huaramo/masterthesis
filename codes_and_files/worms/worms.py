from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import math
import itertools
import argparse
import collections
import warnings
import re
import string
warnings.filterwarnings("ignore")

def read_fasta(file_name):
    '''
    Read and parse FASTA-file.
    :param file_name: The name of FASTA-file.
    '''
    #yield is an efficient alternative of return statement.
    for record in SeqIO.parse(file_name, 'fasta'):
        yield [str(record.id).split('|')[1], str(record.seq).replace(' ', '')]

def seq_identity(seq1, seq2):
    '''
    Compute the sequence identity of two sequences.
    :param seq1: sequence, of course.
    '''
    #Replace the characters from U to X because there is no such characters in BLOSUM62.
    seqs=[seq1.replace('U', 'X'), seq2.replace('U', 'X')]
    
    #Global pairwise alignment with BLOSUM62, gap opening penalty 11 and gap extension penalty 1.
    align= pairwise2.align.globaldd(seqs[0], seqs[1], matlist.blosum62, -11, -1, -11, -1)

    #There are many possible results of the alignment. Choose the one with the best score.
    align_strs=pd.DataFrame(align).sort_values(2, ascending=False).iloc[0].values[0:2]
    
    #Pack the aligned sequences into a DataFrame of characters.
    df_seqs=pd.DataFrame(np.asarray([list(align_strs[0]), list(align_strs[1])]).T.tolist())

    #The number of the matches can thus be calculated efficiently without For-loop.
    match = len(df_seqs.loc[df_seqs[0]==df_seqs[1]])

    #The length of the longer sequence, then implement the formula in return statement.
    max_pos=np.argmax([len(seqs[0]), len(seqs[1])])
    return match/len(seqs[max_pos])

                      
def add_seqidens(df):
    '''
    Compute sequence identity and add this column to the temporary DataFrame.
    :param df: Main DataFrame that should be updated with the column of sequence identity.
    '''
    
    df['seq_iden']=df.apply(lambda x: seq_identity(x['seq_x'], x['seq_y']), axis=1)
    return df

def add_corrs(df, cols1, cols2):
    '''
    Compute expression correlations and add this column to the temporary DataFrame.
    :param df: The DataFrame to be updated.
    :param cols1: The expression profile of the gene 1. cols2 has the same meaning.
    '''
    #Compute Pearson correlation coefficient and add this column to the temporary DataFrame.
    df['p_corr']=df.apply(lambda x: stats.pearsonr(x[cols1], x[cols2])[0], axis=1)

    #The following lines are the rank-score similarity proposed by Chen et al.(2012).
    #Remove different suffixes of all columns, preserve the name of the organ, so that the subtraction can be performed at the end.
    cols=list(map(lambda x: x.split('_')[1], cols1))
    r1=df[cols1].rank(axis=0,pct=True)
    r1.columns=cols
    r2=df[cols2].rank(axis=0,pct=True)
    r2.columns=cols
    
    #Compute rank-score similarity, then add the column into temporary DataFrame.
    df['r_sim']=1-(r1.div(r1.abs().sum(axis=1), axis=0)-r2.div(r2.abs().sum(axis=1), axis=0)).abs().sum(axis=1)
     
    return df        


def parse_homology(homolog_file, homolog_type):
    '''
    Parse the homolgo file from OrthoVenn2 into DataFrame.
    :param homolog_file: The homolog file from OrthoVenn2.
    :param homolog_type: The type of homolog assigned as string.
    '''
    #Read in the homolog file and drop the unnecessary columns.
    homolog=pd.read_csv(homolog_file, sep='\t', header=None).drop(2, axis=1)

    #Change the column names and assign the corresponding homolog file.
    homolog.columns=['p1', 'p2']
    homolog['type']=homolog_type

    #Extract only the Ensembl IDs in every columns.
    homolog['p1']=homolog['p1'].apply(lambda x: x.split('|')[1])
    homolog['p2']=homolog['p2'].apply(lambda x: x.split('|')[1])
    return homolog

def merge_homology(df, origin_col, change_col, homology):
    '''
    Merge the DataFrame to the homolog file.
    :param df: DataFrame to be merged to homolog file.
    :param origin_col: First column name of the homolog pair.
    :param change_col: second column name of the homolog pair.
    :param homology: The homolog file.
    '''
    tmp_df=df.merge(homology, on=origin_col)
    cols_selected=[change_col]
    for c in df.columns[1:]:
        cols_selected.append(c)
    df.columns=cols_selected

    #The rows with the same gene and duplicated rows are droped in the following lines.  
    tmp_df=df.merge(tmp_df, on=change_col).drop_duplicates()
    select=tmp_df[origin_col]!=tmp_df[change_col]
    tmp_df=tmp_df.loc[select]
    return tmp_df

def select_rnaseq(rnaseq):
    '''
    Filter the genes with unreasonably high mean and low variance.
    :param rnaseq: RNA-seq dataset.
    '''
    rnaseq['mean']=rnaseq.mean(axis=1).values 
    select=(rnaseq['mean'] <= 9000) 
    rnaseq=rnaseq.loc[select] 
    rnaseq['var']=rnaseq.var(axis=1).values 
    select=(rnaseq['var'] >= 2) 
    return select

def main():
    '''
    The dataset of worms is more complicate to handle for the following reasons because
    it only uses the transcript IDs of the C.elegans for representing the genes of all of the worms.
    
    '''


    #Read in the RNA-seq dataset.
    rnaseq=pd.read_csv('./worms/rnaseq_worms.tsv', sep='\t', index_col=0).astype(float).reset_index()
    
    #Parse the homolog file and assign the type of homolog.
    para=parse_homology('./worms/ortho_inparalogs.txt', 'paralog')
    orth=parse_homology('./worms/ortho_orthologs.txt', 'ortholog') 
    co_orth=parse_homology('./worms/ortho_co_orthologs.txt', 'ortholog')

    #Read in the two mapping files and draw the relevant columns.
    mapping=pd.read_csv('./worms/ele_uniprot.tsv', sep='\t')[['gene_stable_id', 'transcript_stable_id']]
    #Change the column names
    mapping.columns=['p2', 'GENEID']

    #Remove the variety of suffixes that cause confusion in the mapping file.
    for alphabet in list(string.ascii_lowercase):
        mapping['GENEID']=mapping['GENEID'].apply(lambda x: x.split(alphabet)[0])
    mapping['GENEID']=mapping['GENEID'].apply(lambda x: '.'.join(x.split('.')[0:2]))

    #Merge the mapping file to the RNA-seq dataset.
    rnaseq=mapping.merge(rnaseq, on ='GENEID').drop('GENEID',axis=1)

    #Merge the ortholog file to the RNA-seq dataset. This can retrieve the corresponding C. briggsae IDs. 
    rnaseq=orth.merge(rnaseq, on='p2').drop_duplicates()

    #Edit the column names of the DataFrame, so that the mean of gene expression values of every developmental stage can be calculated. 
    cols=['p1', 'p2', 'type']  
    main_cols=list(map(lambda x: '_'.join(x.split('_')[0:2]), rnaseq.columns[3:].tolist()))
    cols=np.concatenate([cols, main_cols])
    rnaseq.columns=cols
    cols_selected=['cbr_E', 'cbr_L1', 'cbr_ref', 'cel_E', 'cel_L1', 'cel_ref']
    for c in cols_selected:
        rnaseq[c+'_mean']=rnaseq[c].mean(axis=1)

    #Drop the original selected columns of the original DataFrame.
    rnaseq=rnaseq.drop(cols_selected, axis=1)

    #Separate the columns of species C. elegans and C. briggsae.
    cols_cel=list(filter(lambda x: x[:3] == 'cel', rnaseq.columns))
    cols_cbr=list(filter(lambda x: x[:3] == 'cbr', rnaseq.columns))

    #Filter the genes with unreasonably high mean and low variance.
    select_cel=select_rnaseq(rnaseq[cols_cel])
    select_cbr=select_rnaseq(rnaseq[cols_cbr])
    select=select_cel & select_cbr
    rnaseq=rnaseq.loc[select].drop_duplicates().dropna()
    rnaseq_orth=rnaseq.copy() 
    
    #Drop the information of the pair partners. 
    rnaseq_cbr=rnaseq.drop(['p2', 'type'], axis=1).drop(cols_cel, axis=1)
    para_cbr=merge_homology(rnaseq_cbr, 'p1', 'p2', para)
    rnaseq_cel=rnaseq.drop(['p1', 'type'], axis=1).drop(cols_cbr, axis=1)
    para_cel=merge_homology(rnaseq_cel, 'p2', 'p1', para)
    
    #Preparing for the next merging process for co-orthologs and paralogs.
    cols=['p2']
    for c in rnaseq_cbr.columns[1:]:
        if c[:3] == 'cbr':
            cols.append(c.split('_')[1])
        else:
            cols.append(c)
    rnaseq_cbr.columns=cols
    rnaseq_cel.columns=cols
    rnaseq=pd.concat([rnaseq_cbr, rnaseq_cel])
    
    #Preparing the two DataFrame for analyses: ortholog RNA-seq and paralog-co-ortholog RNA-seq DataFrames.
    #But please note that co-orthologs are counted as orthologs in the later analyses. 
    co_orth_rnaseq=merge_homology(rnaseq, 'p2', 'p1', co_orth)
    para_cbr.columns=co_orth_rnaseq.columns
    para_cel.columns=co_orth_rnaseq.columns
    rnaseq=pd.concat([co_orth_rnaseq, para_cbr, para_cel])
    prot=pd.DataFrame(read_fasta('./worms/all_protein_seq.fasta'), columns=['p2', 'seq'])
    rnaseq_orth=merge_homology(prot, 'p2', 'p1', rnaseq_orth).drop_duplicates().dropna() 
    rnaseq=merge_homology(prot, 'p1', 'p2', rnaseq).drop_duplicates().dropna()
    
    #Sample the same number of orthologs and paralogs, then compute the correlations and sequence identity for every portion of the DataFrame.
    rnaseq_orth1=rnaseq_orth.sample(n=1916)
    rnaseq_orth1=add_corrs(rnaseq_orth1, cols_cel, cols_cbr)
    rnaseq_orth1=add_seqidens(rnaseq_orth1)
    cols1=list(filter(lambda x: x.split('_')[-1]=='x', rnaseq.columns))
    cols2=list(filter(lambda x: x.split('_')[-1]=='y', rnaseq.columns))
    rnaseq_orth2=rnaseq.loc[rnaseq['type']=='ortholog']
    rnaseq_orth2=add_corrs(rnaseq_orth2, cols1[1:], cols2[1:])
    rnaseq_orth2=add_seqidens(rnaseq_orth2)
    rnaseq_para=rnaseq.loc[rnaseq['type']=='paralog']
    rnaseq_para=add_corrs(rnaseq_para, cols1[1:], cols2[1:])
    rnaseq_para=add_seqidens(rnaseq_para)

    #Select the columns that need to be written into the TSV-file
    select=['type','p_corr', 'r_sim', 'seq_iden']
    rnaseq=pd.concat([rnaseq_orth1[select], rnaseq_orth2[select], rnaseq_para[select]])
    rnaseq.to_csv('./worms/summary_be.tsv', sep='\t')
           
main()


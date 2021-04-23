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
    '''Read and parse FASTA-file.
    :param file_name: The name of FASTA-FILE.
    '''
    for record in SeqIO.parse(file_name, 'fasta'):
        yield [str(record.id).split('|')[1], str(record.seq).replace(' ', '')]

               
def seq_identity(seq1, seq2):
    '''
    Compute the sequence identity of the two sequences.
    :param seq1: sequence, of course.
    '''
    #Replace the characters from U to X because there is no such characters in BLOSUM62.
    seqs=[seq1.replace('U', 'X'), seq2.replace('U', 'X')]

    #Global pairwise alignment with BLOSUM62, gap opening penalty 11 and gap extension penalty 1.
    align= pairwise2.align.globaldd(seqs[0], seqs[1], matlist.blosum62, -11, -1, -11, -1)

    #There are many possible results of the alignment. Choose the one with the best score.
    align_strs=pd.DataFrame(align).sort_values(2, ascending=False).iloc[0].values[0:2]
    
    #Pack the aligned sequences into a DataFrame.
    df_seqs=pd.DataFrame(np.asarray([list(align_strs[0]), list(align_strs[1])]).T.tolist())

    #The number of matches can thus be calculated efficiently without For-loop.
    match = len(df_seqs.loc[df_seqs[0]==df_seqs[1]])

    #The length of the longer sequence. 
    max_pos=np.argmax([len(seqs[0]), len(seqs[1])])

    #Implement the formula in return statement.
    return match/len(seqs[max_pos])

def add_seqidens(df):
    '''
    Compute sequence similarity and add this column to the temporary DataFrame.
    :param df: Main DataFrame that needs to be updated with the column of sequence identity.
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
    #Ranked the expression values column-wise.
    r1=df[cols1].rank(axis=0,pct=True)
    r1.columns=cols
    r2=df[cols2].rank(axis=0,pct=True)
    r2.columns=cols

    #Implement the formula of the rank-score similarity.
    df['r_sim']=1-(r1.div(r1.abs().sum(axis=1), axis=0)-r2.div(r2.abs().sum(axis=1), axis=0)).abs().sum(axis=1)
     
    return df        
  
def parse_homology(homolog_file, homolog_type):
    '''
    Parse the homolog file from OrthoVenn2 into DataFrame.
    :param homolog-file: The homolog file from OrthoVenn2.
    :param homolog type: The type of homolog assigned as string.
    '''
    #Read in the homolog file.
    homolog=pd.read_csv(homolog_file, sep='\t', header=None).drop(2, axis=1)

    #Change the column names and assign the corresponding homolog type.
    homolog.columns=['p1', 'p2']
    homolog['type']=homolog_type

    #Extract only the Ensembl IDs in every column.
    homolog['p1']=homolog['p1'].apply(lambda x: x.split('|')[1])
    homolog['p2']=homolog['p2'].apply(lambda x: x.split('|')[1])
    return homolog

def merge_homology(df, origin_col, change_col, homology):
    '''
    Merge the DataFrame to the homolog file.
    :param df: DataFrame to be merged to homolog file.
    :param origin_col: First column name of the homolog pair.
    :param change_col: Second column name of the homolog pair.
    :param homology: The homolog file.
    '''
    tmp_df=df.merge(homology, on=origin_col)
    cols_selected=[change_col]
    for c in df.columns[1:]:
        cols_selected.append(c)
    df.columns=cols_selected
    tmp_df=df.merge(tmp_df, on=change_col).drop_duplicates()
    select=tmp_df[origin_col]!=tmp_df[change_col]
    tmp_df=tmp_df.loc[select]
    return tmp_df


def main():
    #Read in the RNA-seq dataset.
    rnaseq=pd.read_csv('./plants/rnaseq_plants.tsv', sep='\t', index_col=0)

    #Drop the undesired expression vectors, either variance = 0 or unreasonable high mean.
    rnaseq['mean']=rnaseq.mean(axis=1).values
    select=(rnaseq['mean'] <= 9000)
    rnaseq=rnaseq.loc[select]
    rnaseq['var']=rnaseq.var(axis=1).values
    select=(rnaseq['var'] >= 2) 

    #Drop the help columns of mean and variance.
    rnaseq=rnaseq.loc[select].drop(['mean', 'var'], axis=1).reset_index()

    #Change the column name of the IDs to p2.
    cols=['p2']
    for c in rnaseq.columns[1:]:
        cols.append(c)
    rnaseq.columns=cols
    
    #Read in the protein FASTA-file and pack into a DataFrame.
    prot=pd.DataFrame(read_fasta('./plants/all_protein_seq.fasta'), columns=['p2', 'seq'])
    rnaseq=prot.merge(rnaseq, on='p2')
    
    #Parse the homolog file, assign the type of the homolog and pack them into a DataFrame.
    para=parse_homology('./plants/ortho_inparalogs.txt', 'paralog')
    orth=parse_homology('./plants/ortho_orthologs.txt', 'ortholog')
    co_orth=parse_homology('./plants/ortho_co_orthologs.txt', 'ortholog')
    homologs=pd.concat([orth, co_orth, para])
    rnaseq=merge_homology(rnaseq, 'p2', 'p1', homologs)
    
    #The following lines draw the same number of orthologs and paralogs. 
    orth=rnaseq.loc[(rnaseq['type']=='ortholog')]
    para=rnaseq.loc[(rnaseq['type']=='paralog')]
    if len(orth.index) <= len(para.index):
        length=len(orth.index) 
    else: 
        length=len(para.index) 
    orth=orth.sample(n=length) 
    para=para.sample(n=length)
    rnaseq=pd.concat([orth, para])

    #Separate the columns of every homolog pair for correlation analysis.
    cols1=list(filter(lambda x: x.split('_')[-1]=='x', rnaseq.columns))
    cols2=list(filter(lambda x: x.split('_')[-1]=='y', rnaseq.columns))
    
    #Compute the expression correlations and add columns to temporary DataFrame.  
    rnaseq=add_corrs(rnaseq, cols1[1:], cols2[1:])
    rnaseq=add_seqidens(rnaseq)

    #Remove the NaNs and write the TSV-file.
    rnaseq.dropna().to_csv('./plants/summary_ag.tsv', sep='\t') 
           
main()


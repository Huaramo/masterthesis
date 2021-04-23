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
    Read and parse the FASTA-file.
    :param file_name: The name of FASTA-file.
    '''
    #yield is an efficient alternative of return statement.
    for record in SeqIO.parse(file_name, 'fasta'):
        yield [str(record.id).split('|')[1], str(record.seq).replace(' ', '')]

                
def seq_identity(seq1, seq2):
    '''
    Compute the sequence identity for a pair of sequences.
    :param seq1: sequence, of course. 
    '''
    #Replace the characters from * and U to X because there are no such characters in BLOSUM62 matrix.
    seqs=[seq1.replace('*', 'X').replace('U', 'X'), seq2.replace('*', 'X').replace('U', 'X')]

    #Pairwise global alignment with BLOSUM62 matrix, gap opening penalty 11 and gap extension penalty 1.
    align= pairwise2.align.globaldd(seqs[0], seqs[1], matlist.blosum62, -11, -1, -11, -1)
    #There are many possible results of the alignments. Choose the one with the highest score.
    align_strs=pd.DataFrame(align).sort_values(2, ascending=False).iloc[0].values[0:2]
    #Pack the aligned sequences into a DataFrame of characters.
    df_seqs=pd.DataFrame(np.asarray([list(align_strs[0]), list(align_strs[1])]).T.tolist())
    #The number of matches can be calculated efficiently without For-loop.
    match = len(df_seqs.loc[df_seqs[0]==df_seqs[1]])
    #The length of the longer sequence. Implement the formula in return statement.
    max_pos=np.argmax([len(seqs[0]), len(seqs[1])])
    
    return match/len(seqs[max_pos])



def add_seqidens(df):
    '''
    Compute sequence identity and add this column to the temporary DataFrame.
    :param df: Main DataFrame that needs to be updated with the column of the sequence identity.
    '''
    
    df['seq_iden']=df.apply(lambda x: seq_identity(str(x['seq_x']), str(x['seq_y'])), axis=1)
    return df

def add_corrs(df, cols1, cols2):
    '''Update the main DataFrame with expression correlation coefficients.'''
    df['p_corr']=df.apply(lambda x: stats.pearsonr(x[cols1], x[cols2])[0], axis=1)
    cols=list(map(lambda x: x.split('_')[1], cols1))
    r1=df[cols1].rank(axis=0,pct=True)
    r1.columns=cols
    r2=df[cols2].rank(axis=0,pct=True)
    r2.columns=cols
    df['r_sim']=1-(r1.div(r1.abs().sum(axis=1), axis=0)-r2.div(r2.abs().sum(axis=1), axis=0)).abs().sum(axis=1)
     
    return df        


   
    
def parse_homology(homolog_file, homolog_type):
    '''
    Parse the homolog file from OrthoVenn2 into DataFrame.
    :param homolog_file: The homolog file from OrthoVenn2.
    :param homolog_type: The type of homolog assigned as string.
    '''
    #Read in the homolog file and drop the unnecessary columns.
    homolog=pd.read_csv(homolog_file, sep='\t', header=None).drop(2, axis=1)
    #Change the column names and assign the corresponding homolog type.
    homolog.columns=['p1', 'p2']
    homolog['type']=homolog_type
    #Extract only the Ensembl IDs in every column:
    homolog['p1']=homolog['p1'].apply(lambda x: x.split('|')[1])
    homolog['p2']=homolog['p2'].apply(lambda x: x.split('|')[1])
    return homolog

def merge_homology(df, origin_col, change_col, homology):
    '''
    Merge the DataFrame to the homolog file
    :param df: DataFrame to be merged to homolog file.
    :param origin_col: The first column name of the homolog pair.
    :param change_col: The second column name of the homolog pair.
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


def parsing():
    '''Initialize parser object.'''
    parser = argparse.ArgumentParser()
    '''Add argument'''
    parser.add_argument('--fasta', type=str, help='Input paths of FASTA-file with 2 species.')
    parser.add_argument('--rnaseq', type=str, help='Input paths of preprocessed RNA-seq with 2 species.')
    parser.add_argument('--homologs', type=str, nargs=3, help='Input paths of inparalogs, orthologs, co-orthologs defined by the OrthoVenn2.')
    parser.add_argument('--out', type=str, help='Output path of the summary file.')
    '''Parse the arguments, then we can access every input parameter by using the name of the argument.'''
    args = parser.parse_args()
    return args

def main():
    args=parsing()
    rnaseq=pd.read_csv(args.rnaseq, sep='\t', index_col=0)###
    
    #Drop the undesired expression vectors, either variance = 0 or unreasonable high mean. 
    rnaseq['mean']=rnaseq.mean(axis=1).values
    select=(rnaseq['mean'] <= 9000)
    rnaseq=rnaseq.loc[select]
    rnaseq['var']=rnaseq.var(axis=1).values
    select=(rnaseq['var'] >= 2)
    #Drop the help columns.
    rnaseq=rnaseq.loc[select].drop(['mean', 'var'], axis=1).reset_index()
    
    
    cols=['p2']
    for c in rnaseq.columns[1:]:
        cols.append(c)
    rnaseq.columns=cols
    
    prot=pd.DataFrame(read_fasta(args.fasta), columns=['p2', 'seq'])
    rnaseq=prot.merge(rnaseq, on='p2')
    para=parse_homology(args.homologs[0], 'paralog')
    orth=parse_homology(args.homologs[1], 'ortholog')
    co_orth=parse_homology(args.homologs[2], 'ortholog')
    homologs=pd.concat([orth, co_orth, para])
    rnaseq=merge_homology(rnaseq, 'p2', 'p1', homologs)

    #The following lines draw the same amount of orthologs and paralogs for fair comparison.
    orth=rnaseq.loc[(rnaseq['type']=='ortholog')]
    para=rnaseq.loc[(rnaseq['type']=='paralog')]
    
    if len(orth.index) <= len(para.index):
        length=len(orth.index) 
    else: 
        length=len(para.index) 
    orth=orth.sample(n=length) 
    para=para.sample(n=length)     
    #print(len(orth.index))
    #print(len(para.index))
    rnaseq=pd.concat([orth, para])
    
    #Keep the necessary columns for correlation analysis        
    cols1=list(filter(lambda x: x.split('_')[-1]=='x', rnaseq.columns))
    cols2=list(filter(lambda x: x.split('_')[-1]=='y', rnaseq.columns))

    #Compute expression correlation and sequence identity.
    rnaseq=add_corrs(rnaseq, cols1[1:], cols2[1:])
    rnaseq=add_seqidens(rnaseq)
    rnaseq.dropna().to_csv(args.out, sep='\t')
    
           
main()


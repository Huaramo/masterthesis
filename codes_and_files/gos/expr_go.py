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

def resnik(com_ancs, probs):
    '''
    Compute the Resnik semantic similarity.
    :param com_ancs: The set of common ancestor GO-terms.
    :param probs: The table of conditional probability and relative frequency for individual ontology.
    '''
    #Select the part of probability table that has the information of common ancestors.
    select = pd.Index(com_ancs).intersection(probs.index)
    probs=probs.loc[select]

    #If there is no information of common ancestors, give resnik similarity the value 0.    
    if len(probs) == 0:
        return 0
    #Implement the formula directly.        
    return np.amax(-np.log2(probs['rel_freq']))

                
def sum_info_gain(probs, arr):
    '''
    Compute the sum of information accretions.
    :param probs: The table of conditional probability and relative frequency for individual ontology.
    :param arr: This array is a set of propagated annotations(GO-terms).
    '''
    #Pack the set of propagated annotations into a DataFrame with only one column.
    df_arr=pd.DataFrame(arr, index=np.arange(len(arr)), columns=['node'])

    #Take the relevant part of the probability table.
    df_arr=df_arr.merge(probs, on='node')
    summation=np.sum(-np.log2(df_arr['cond_prob']))
    
    return summation

def lin(row):
    '''
    Compute Lin semantic similarity, a help function in schlicker().
    :param row: The row of probability table.()
    '''
    #Implement directly the formula.
    return -2*row['resnik']/(np.log2(row['rel_freq_x'])+np.log2(row['rel_freq_y']))

def maryland_bridge(gos1, gos2):
    '''
    Compute Maryland-Bridge functional similarity.
    :param gos1: The set of propagated annotations. gos2 has the same meaning.
    '''
    #Split the propagated annotations which are connected by , .
    gos1=gos1.split(',')
    gos2=gos2.split(',')

    #There are somehow nan, which must be filtered, in the string..
    gos1=list(filter(lambda p: p != 'nan', gos1))
    gos2=list(filter(lambda p: p != 'nan', gos2))

    #Intersection of the two sets of propagated annotations.
    anb=np.intersect1d(gos1, gos2, assume_unique=True)

    #Implement directly the formula. 
    return (len(anb)/(2*len(gos1)))+(len(anb)/(2*len(gos2)))


def schlicker(gos1, gos2, probs):
    '''
    Compute Schlicker functional similarity.
    :param gos1: The set of gene annotations. gos2 has the same meaning.
    :param probs: The probability table of conditional probability and relative frequency for individual ontology.
    '''
    #Split the propagated annotations which are connected by , .    
    gos1=gos1.split(',') 
    gos2=gos2.split(',') 

    #Filter out the NaNs.    
    select1=list(filter(lambda p: p != 'nan', gos1))
    select2=list(filter(lambda p: p != 'nan', gos2))
    probs=probs.set_index('node')

    #The long if-condition is to prevent the case of empty table after the selection process.
    if len(probs.index.intersection(pd.Index(select1))) > 0 and len(probs.index.intersection(pd.Index(select2))) > 0:
        #The relevant part of probability table is selected.
        probs_selected1=probs.loc[select1].dropna().reset_index()
        probs_selected2=probs.loc[select2].dropna().reset_index()
        
        #Perform the cross-join. 
        probs_selected1['key']=1
        probs_selected2['key']=1
        probs_selected=probs_selected1.merge(probs_selected2, on='key')

        #Compute and add the column of common ancestors on the selected table. 
        probs_selected['com_ancs']=probs_selected[['path_x', 'path_y']].apply(lambda x: np.intersect1d(x['path_x'].split(','), x['path_y'].split(','), assume_unique=True).tolist(), axis=1)

        #Compute Resnik and Lin semantic similarity and add these columns to the selected DataFrame.
        probs_selected['resnik']=probs_selected['com_ancs'].apply(lambda x: resnik(x, probs)) 
        probs_selected['lin']=probs_selected.apply(lambda x: lin(x), axis=1)

        #Find the best matches of every GO-term from both genes.
        tmp1=probs_selected.groupby('node_x')['lin'].max().values.tolist()
        tmp2=probs_selected.groupby('node_y')['lin'].max().values.tolist()

        #Implement the formula.
        return (np.sum(tmp1)+np.sum(tmp2))/(len(gos1)+len(gos2))
    else:
        #If the selection yields empty table, then pass NaN.
        return float('nan')

 
def yang_clark(gos1, gos2, probs):
    '''
    Compute Schlicker functional similarity.
    :param gos1: The set of gene annotations. gos2 has the same meaning.
    :param probs: The probability table of conditional probability and relative frequency for individual ontology.
    '''

    #Split the propagated annotations which are connected by , .
    gos1=gos1.split(',')
    gos2=gos2.split(',')
    
    #Filter NaN after the split.
    gos1=list(filter(lambda p: p != 'nan', gos1))
    gos2=list(filter(lambda p: p != 'nan', gos2))

    #Compute the two set differences and union.
    #anb=np.intersect1d(gos1, gos2, assume_unique=True)
    a_b = np.setdiff1d(gos1, gos2, assume_unique=True)
    b_a = np.setdiff1d(gos2, gos1, assume_unique=True)
    aub = np.union1d(gos1, gos2)
    denom = sum_info_gain(probs, aub)
    #This case is to prevent division by zero.
    if denom != 0:
        #Implement the formula.                                                                                                                                                                 #print('enter4')
        yc = 1 - (math.sqrt(sum_info_gain(probs, a_b) ** 2 + sum_info_gain(probs, b_a) ** 2) / sum_info_gain(probs, aub))
    else:
        #Assign NaN.
        yc = float('nan')
    
    
    return yc



def add_semsims(df_tmp, anc):
    '''
    Compute semantic similarity and add this column to the temporary DataFrame.
    :param df_tmp: Main DataFrame that should be updated with the column of sequence identity
    :param anc: Root node of the ontology, BP, CC or MF.
    '''
    #Read in the file then drop NaNs.
    probs=pd.read_csv('./gos/probs_'+anc+'.tsv', sep='\t', index_col=0).dropna()        

    #Import the paths for all GO-terms      
    paths=pd.read_csv('./gos/paths_prob_'+anc+'.tsv', sep='\t', index_col=0).dropna()
    
    #Change the column name and merge from both sides to construct the table with ortho-/paralog pairs
    cols=['xref_y']
    for c in paths.columns[1:]:
        cols.append(c)
    paths.columns=cols
    df_tmp=paths.merge(df_tmp, on='xref_y')
    cols=['xref_x']
    for c in paths.columns[1:]:
        cols.append(c)
    paths.columns=cols
    df_tmp=paths.merge(df_tmp, on='xref_x').drop_duplicates().dropna()

    #Draw the same number of ortho-/paralogs for fair comparison.
    tmp_orth = df_tmp.loc[(df_tmp['type']=='ortholog')]
    tmp_para = df_tmp.loc[(df_tmp['type']=='paralog')]
    if len(tmp_orth.index) <= len(tmp_para.index):
        length=len(tmp_orth.index)
    else:
        length=len(tmp_para.index)

    tmp_orth=tmp_orth.sample(n=length)
    tmp_para=tmp_para.sample(n=length)
    df_tmp=pd.concat([tmp_orth, tmp_para]).reset_index().drop('index', axis=1)
    
    #Compute Schlicker similarity
    df_tmp['sc_sim']=df_tmp.apply(lambda x: schlicker(x['go_term_x'], x['go_term_y'], probs), axis=1)
    #Compute Yang-Clark similarity
    df_tmp['yc_sim']=df_tmp.apply(lambda x: yang_clark(x['path_x'], x['path_y'], probs), axis=1)
    #Compute Maryland-bridge similarity
    #df_tmp['mb_sim']=df_tmp.apply(lambda x: maryland_bridge(x['path_x'], x['path_y']), axis=1)

    return df_tmp

def add_corrs(df, cols1, cols2):
    '''
    Compute expression correlations and add this column to the temporary DataFrame.
    :param df: The DataFrame to be updated.
    :param cols1: The expression profile of the gene 1. cols2 has the same meaning.  
    '''
    #Compute the expression correlation 
    df['p_corr']=df.apply(lambda x: stats.pearsonr(x[cols1], x[cols2])[0], axis=1)
    #df['sp_corr']=df.apply(lambda x: stats.spearmanr(x[cols1], x[cols2])[0], axis=1)
    #df['tau_corr']=df.apply(lambda x: stats.kendalltau(x[cols1], x[cols2])[0], axis=1)
    cols=list(map(lambda x: x.split('_')[1], cols1))

    #The following lines compute rank-score similarity proposed by Chen et al.(2012).
    #Step 1: Rank the expression values columnwise and convert the ranks into percentiles.
    r1=df[cols1].rank(axis=0,pct=True)
    r1.columns=cols
    r2=df[cols2].rank(axis=0,pct=True)
    r2.columns=cols

    #Compute the rank-score similarity.
    df['r_sim']=1-(r1.div(r1.abs().sum(axis=1), axis=0)-r2.div(r2.abs().sum(axis=1), axis=0)).abs().sum(axis=1)
     
    return df        
    
def parse_homology(homolog_file, homolog_type):
    '''
    #Parse the homolog file from OrthoVenn2 into DataFrame.
    :param homolog_file: The homolog file from OrthoVenn2
    :param homolog_type: The type of homolog assigned as string.
    '''
    #Read in the homolog file and drop the unnecessary column.
    homolog=pd.read_csv(homolog_file, sep='\t', header=None).drop(2, axis=1)
    homolog.columns=['p1', 'p2']

    #Assign the corresponding homolog type.
    homolog['type']=homolog_type
    
    #Extract only the Ensembl-IDs.
    homolog['p1']=homolog['p1'].apply(lambda x: x.split('|')[1])
    homolog['p2']=homolog['p2'].apply(lambda x: x.split('|')[1])
    return homolog

def merge_homology(df, origin_col, change_col, homology):
    '''
    Merge the DataFrame to the homology file.
    :param df: DataFrame to be merged to homolog file.
    :param origin_col: First column name of the homolog pair.
    :param change_col: Second column name of the homolog pair. 
    :param homology: The homology file. 
    '''
    
    tmp_df=df.merge(homology, on=origin_col)
    cols_selected=[change_col]
    for c in df.columns[1:]:
        cols_selected.append(c)
    df.columns=cols_selected
    #The repetitive rows and the rows in which the both sides have the same information should be dropped.
    tmp_df=df.merge(tmp_df, on=change_col).drop_duplicates()
    select=tmp_df[origin_col]!=tmp_df[change_col]
    tmp_df=tmp_df.loc[select]
    return tmp_df


def parsing():
    '''Initialize parser object.'''
    parser = argparse.ArgumentParser()
    '''Add argument'''
    parser.add_argument('--rnaseq', type=str, help='Input paths of preprocessed RNA-seq with 2 species.')
    parser.add_argument('--maps', type=str, nargs=2, help='Input paths of 2 Uniprot-ID-mapping-files from Ensembl FTP-site.')
    parser.add_argument('--homologs', type=str, nargs=3, help='Input paths of inparalogs, orthologs, co-orthologs defined by the OrthoVenn2.')
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
    rnaseq=rnaseq.loc[select].drop(['mean', 'var'], axis=1).reset_index()
    
    
    cols=['p2']
    for c in rnaseq.columns[1:]:
        cols.append(c)
    rnaseq.columns=cols

    #Parse the homolog files, assign the homolog type and combine them into a DataFrame..
    para=parse_homology(args.homologs[0], 'paralog')
    orth=parse_homology(args.homologs[1], 'ortholog')
    co_orth=parse_homology(args.homologs[2], 'ortholog')
    homologs=pd.concat([orth, co_orth, para])

    #Merge the RNA-seq dataset to the homolog file.
    tmp=merge_homology(rnaseq, 'p2', 'p1', homologs).drop_duplicates().dropna()

    #Read in mapping files.
    map1=pd.read_csv(args.maps[0], sep='\t')
    map2=pd.read_csv(args.maps[1], sep='\t') 
    map=pd.concat([map1, map2]) 

    #Select the relevant columns: Uniprot-IDs and Ensembl-IDs. Then merge these columns to the temporary DataFrame. 
    map=map[['xref', 'gene_stable_id']] 
    map.columns=['xref', 'p2'] 
    tmp=map.merge(tmp, on='p2') 
    map.columns=['xref', 'p1'] 
    tmp=map.merge(tmp, on='p1')
    cols1=list(filter(lambda x: x.split('_')[-1]=='x', tmp.columns))
    cols2=list(filter(lambda x: x.split('_')[-1]=='y', tmp.columns))
    
    #For every ontology, compute the functional simiarity and expression correlation measures and add columns to the temporary DataFrame.
    for anc in ['bp', 'cc', 'mf']:
        tmp_copy=tmp.copy()
        tmp_copy=add_semsims(tmp_copy, anc)        
        tmp_copy=add_corrs(tmp_copy, cols1[2:], cols2[2:])
        tmp_copy.dropna().to_csv('./gos/summary_expr_go_'+args.rnaseq.split('.')[1].split('/')[2]+'_'+anc+'.tsv', sep='\t')
main()


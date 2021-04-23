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

    #There are somehow nan, which must be filtered, in the string... 
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
    
    #Prepare the selection on table.
    select1=list(filter(lambda p: p != 'nan', gos1))
    select2=list(filter(lambda p: p != 'nan', gos2))
    probs=probs.set_index('node')

    #The long if-condition is to prevent the case of empty table after the selection process.
    if len(probs.index.intersection(pd.Index(select1))) > 0 and len(probs.index.intersection(pd.Index(select2))) > 0:
        #The relevant part of probability table is selected. 
        probs_selected1=probs.loc[select1].dropna().reset_index()
        probs_selected2=probs.loc[select2].dropna().reset_index()

        #Perform cross-join on the selected table.
        probs_selected1['key']=1
        probs_selected2['key']=1
        probs_selected=probs_selected1.merge(probs_selected2, on='key')

        #Compute and add the column of common ancestors on the selected table.
        probs_selected['com_ancs']=probs_selected[['path_x', 'path_y']].apply(lambda x: np.intersect1d(x['path_x'].split(','), x['path_y'].split(','), assume_unique=True).tolist(), axis=1) 

        #Compute and add the column of Resnik semantic similarity.
        probs_selected['resnik']=probs_selected['com_ancs'].apply(lambda x: resnik(x, probs))

        #Compute and add the column of Lin semantic similarity.
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
    a_b = np.setdiff1d(gos1, gos2, assume_unique=True)
    b_a = np.setdiff1d(gos2, gos1, assume_unique=True)
    aub = np.union1d(gos1, gos2)

    #Compute the denominator of yang-clark formula.
    denom = sum_info_gain(probs, aub)

    #Prevent the case of division by zero.
    if denom != 0:
        #Implement the formula.       
        yc = 1 - (math.sqrt(sum_info_gain(probs, a_b) ** 2 + sum_info_gain(probs, b_a) ** 2) / denom)
    else:
        yc = float('nan')
    
    return yc

def read_fasta(file_name):
    '''
    Read and parse FASTA-file.
    :param file_name: The name of the FASTA-file.
    '''
    #yield is an efficient alternative of return statement.
    for record in SeqIO.parse(file_name, 'fasta'):
        yield [str(record.id).split('|')[1], str(record.seq).replace(' ', '')]

def seq_identity(seq1, seq2):
    '''
    Compute the sequence identity of two sequences.
    :param seq1: sequence, of course. 
    '''
    #Replace the characters from * and U to X because there is no such characters in BLOSUM62. 
    seqs=[seq1.replace('*', 'X').replace('U', 'X'), seq2.replace('*', 'X').replace('U', 'X')]
    #Global pairwise alignment with BLOSUM62, gap opening penalty 11 and gap extension penalty 1. 
    align= pairwise2.align.globaldd(seqs[0], seqs[1], matlist.blosum62, -11, -1, -11, -1)
    #There are many possible results of the alignment. Choose the one with the best score. 
    align_strs=pd.DataFrame(align).sort_values(2, ascending=False).iloc[0].values[0:2]
    #Pack the aligned sequences into a DataFrame of characters.
    df_seqs=pd.DataFrame(np.asarray([list(align_strs[0]), list(align_strs[1])]).T.tolist())
    #The number of matches can thus be calculated efficiently without For-loop.
    match = len(df_seqs.loc[df_seqs[0]==df_seqs[1]])
    #The length of the longer sequence then implement the formula in return statement. 
    max_pos=np.argmax([len(seqs[0]), len(seqs[1])])

    return match/len(seqs[max_pos])
       
def add_seqidens(df_tmp, file_fasta):
    '''
    Compute sequence similarity and add this column to the temporary DataFrame.
    :param df_tmp: Main DataFrame that should be updated with the column of sequence identity
    :param file_fasta: File name of the FASTA-file.
    '''
    #Read and parse the FASTA-file and store it into a DataFrame.
    prot=pd.DataFrame(read_fasta(file_fasta))
    prot.columns=['p2', 'seq']
    #Merge the main DataFrame with the FASTA-file.
    df_tmp=prot.merge(df_tmp, on='p2')
    prot.columns=['p1', 'seq']
    #The same as previous merge. 
    df_tmp=prot.merge(df_tmp, on='p1').drop_duplicates().dropna()
    #Compute the sequence identity.
    df_tmp['seq_iden']=df_tmp.apply(lambda x: seq_identity(x['seq_x'], x['seq_y']), axis=1)
    return df_tmp

    

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

    #The following lines draw the same amount of orthologs and paralogs for fair comparison.
    tmp_orth = df_tmp.loc[(df_tmp['type']=='ortholog')]
    tmp_para = df_tmp.loc[(df_tmp['type']=='paralog')]
    if len(tmp_orth.index) <= len(tmp_para.index):
        length=len(tmp_orth.index) 
    else: 
        length=len(tmp_para.index)
    
    tmp_orth=tmp_orth.sample(n=length) 
    tmp_para=tmp_para.sample(n=length)
    #After the extraction, combine them into one DataFrame.
    df_tmp=pd.concat([tmp_orth, tmp_para]).reset_index().drop('index', axis=1)
    
    #Compute Schlicker similarity and add this column to temporary DataFrame.
    df_tmp['sc_sim']=df_tmp.apply(lambda x: schlicker(x['go_term_x'], x['go_term_y'], probs), axis=1)
    #Compute Yang-Clark similarity and add this column to temporary DataFrame.
    df_tmp['yc_sim']=df_tmp.apply(lambda x: yang_clark(x['path_x'], x['path_y'], probs), axis=1)
    #Compute Maryland-bridge similarity and add this column to temporary DataFrame.
    #df_tmp['mb_sim']=df_tmp.apply(lambda x: maryland_bridge(x['path_x'], x['path_y']), axis=1)
        
    return df_tmp

def parse_homology(homolog_file, homolog_type):
    '''
    #Parse the homolog file from OrthoVenn2 into DataFrame.
    :param homolog_file: The homolog file from OrthoVenn2
    :param homolog_type: The type of homolog assigned as string.
    '''
    #Read in the homolog file and drop the unnecessary column.
    homolog=pd.read_csv(homolog_file, sep='\t', header=None).drop(2, axis=1)
    #Change the column names and assign the corresponding homolog type.
    homolog.columns=['p1', 'p2']
    homolog['type']=homolog_type
    #Extract only the Ensembl IDs in every column. 
    homolog['p1']=homolog['p1'].apply(lambda x: x.split('|')[1])
    homolog['p2']=homolog['p2'].apply(lambda x: x.split('|')[1])
    return homolog



def parsing():
    '''Here are the parsers'''

    '''Function Call:
    python yc_sim_hm.py --fastas < 2_FASTA-file> --maps < 2_Uniprot-mapping_file> --go < 1_File_of_GO-term-path_per_protein> --homology < 1_Homology_file>
    --rnaseq <(Optional) 1_File_of_RNA-seq_files> --out <Name_of_the_summary_file>'''


    parser = argparse.ArgumentParser() 
    '''Arguments'''
    parser.add_argument('--fasta', type=str, help='Input paths of FASTA-file with 2 species.')
    parser.add_argument('--maps', type=str, nargs=2, help='Input paths of 2 Uniprot-ID-mapping-files from Ensembl FTP-site.')
    parser.add_argument('--homologs', type=str, nargs=3, help='Input paths of inparalogs, orthologs, co-orthologs defined by the OrthoVenn2.')
     
    
    '''Parse the arguments, then we can access every input parameter by using the name of the argument.'''
    args = parser.parse_args()
    return args

def main():
    args=parsing()
    #Import the homolog files and the Uniprot mapping files and merge them.
    #For every homolog file, assign the type of orthologs
    para=parse_homology(args.homologs[0], 'paralog')
    orth=parse_homology(args.homologs[1], 'ortholog')
    co_orth=parse_homology(args.homologs[2], 'ortholog') 

    #Combine the three homolog files into one large DataFrame.
    homologs=pd.concat([orth, co_orth, para])

    #Read in the ID-mapping files and combine them into one large DataFrame..  
    map1=pd.read_csv(args.maps[0], sep='\t') 
    map2=pd.read_csv(args.maps[1], sep='\t')
    map=pd.concat([map1, map2])

    #Take the necessary columns(Uniprot-IDs and Ensembl-IDs) and merge the table so that it contains the Uniprot-IDs of ortho-/paralog pairs.   
    map=map[['xref', 'gene_stable_id']]
    map.columns=['xref', 'p2'] 
    tmp=homologs.merge(map, on='p2')
    map.columns=['xref', 'p1'] 
    tmp=map.merge(tmp, on='p1')
    #For the 3 ontologies, add the columns of semantic similarity and sequence identity.
    for anc in ['bp', 'cc', 'mf']: 
        tmp_copy=tmp.copy()
        tmp_copy=add_semsims(tmp_copy, anc)
        tmp_copy=add_seqidens(tmp_copy, args.fasta)
        #Write into a file. 
        tmp_copy.to_csv('./gos/summary_gos_'+args.homologs[0].split('_')[-1].split('.')[0]+'_'+anc+'.tsv', sep='\t') 
    
    
main()


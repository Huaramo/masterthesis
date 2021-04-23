#Import all the necessary modules.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import math
import itertools
import argparse
import collections
import warnings
import functools as ft
import time
warnings.filterwarnings("ignore")
               
def selector(item):
    '''
    Select the edge relationships, only the ontological ones.
    It is an help function in child_parent().
    '''
    return item[:4] == 'is_a' or item[:21] == 'relationship: part_of' 


def edit(term):
    '''
    It is an help function in child_parent()
    '''
    edited = []
    for t in term:
        edited.append(t.split(' ! ')[0])

    return edited


def rel_concat(term):
    '''
    It is an help function in child_parent()
    '''
    tmp_str = ''
    for rel in term[2:]:
        tmp_str += rel + ','
    return [term[0], term[1], tmp_str[:-1]]



def real_dfs(graph, start):
    '''
    Perform depth-first traversal of the ontology graph.
    :param graph: The child-parent relationship(ontology graph).
    :param start: Starting GO-term of the depth-first traversal. 
    '''
    #Initialize a recording DataFrame which tells whether a node is visited or not.
    visited = pd.DataFrame(np.zeros(len(graph.index)), index=graph.index)

    #From the starting point, perform depth-first traversal and collect the paths.
    collect = [start]
    if int(visited.loc[start].values[0]) == 0:
        real_dfs1(graph, start, visited, collect)

    return collect


def real_dfs1(graph, ind, visited, collect):
    '''
    Help function of real_df().
    :param graph: Ontology graph.
    :param ind: The GO-term.
    :param visited: The DataFrame recording if a node is visited or not.
    :param collect: A collector of the paths.
    '''
    #Parents of a GO-term.
    neighbors = graph.loc[ind].values[0]

    #If there is no parents for this GO-term, we are at the root node.
    #Give out the path.
    if neighbors == '':
        return collect

    #If not, every parent will be served as a starting point of the DFS. 
    neighbors = neighbors.split(',')

    for n in neighbors:
        if n in visited.index:
            #If a node has not been visited, pack this node into the collector, then mark it as visited. 
            #Keep traversing down until all are visited.
            if int(visited.loc[n].values[0]) == 0:
                collect.append(n)
                visited.at[n] = 1
                real_dfs1(graph, n, visited, collect)


  

def extract_gos(item):
    '''
    Help function in child_parent()
    '''
    elem=item.split(',')
    tmp= []
    for e in elem:
        go=e.split(' ')[-1]
        tmp.append(go)
    tmp=np.unique(tmp)
    tmp_str=''
    for t in tmp:
        tmp_str += t + ','
    return tmp_str[:-1]

def child_parent():
    '''
    Parses the ontology from the OBO-file.
    The ontology is represented as a DataFrame of node and its children. 
    '''
    #Read in and parse the OBO-file. The following lines are too trivial to explain...
    file=open('./gos/go.obo', 'r')
    go_terms=[]
    tmp=[]
    for line in file:
        if line[0:6] == '[Term]':
            go_terms.append(tmp)
            tmp = []
        else:
            tmp.append(line[:-1])
    go_terms=go_terms[1:]
    parent_child_rels=[]
    tmp = []
    for term in go_terms:
        tmp.append(term[0])
        tmp.append(term[2])
        for item in term:
            #selector() chooses the ontological relationships: is_a and part_of.
            if selector(item):
                tmp.append(item)
        parent_child_rels.append(tmp)
        tmp = []
    
    #The following three lines clean up the child-parent relationships as a DataFrame. 
    parent_child_rels=list(map(lambda x: edit(x), parent_child_rels))
    parent_child_rels=list(map(lambda x: rel_concat(x), parent_child_rels))
    parent_child_rels=pd.DataFrame(parent_child_rels).set_index(0)
    
    #The following lines separate the DataFrame into three portions: BP, CC and MF.
    new_index = []
    for ind in parent_child_rels.index:
        new_index.append(ind[4:])
    parent_child_rels.index = new_index
    parent_child_rels = parent_child_rels.replace(np.nan, '', regex=True)
    parent_child_rels[1] = parent_child_rels[1].apply(lambda x: x.split(': ')[1]).replace({'biological_process':'bp', 'cellular_component':'cc', 'molecular_function':'mf'})
    parent_child_rels[2] = parent_child_rels[2].apply(lambda x: extract_gos(x))
    
    #Write the three child-parent DataFrames into TSV-file.   
    ancs=['bp', 'cc', 'mf']
    for anc in ancs:
        parent_child_rels.loc[parent_child_rels[1]==anc].drop(1, axis=1).to_csv('./gos/child_parent_'+anc+'.tsv', sep='\t')
    
    
def paths_go():
    '''
    For every GO-term, find the all paths from itself to its root term, which is called propagated annotations.
    The term all paths simply mean all of the possbile nodes that are passed by DFS and are non-repetitive.
    '''
    
    ancs=['bp', 'cc', 'mf']
    for anc in ancs:
        #Read in the child-parent relationsips 
        child_parent = pd.read_csv('./gos/child_parent_'+anc+'.tsv', sep='\t', index_col=0)
        #Because Pandas interprets the empty string as NaN, this line replace NaN back to empty string for depth-first graph traveral.
        child_parent = child_parent.replace(np.nan, '', regex=True)
        #Write new files which stores all the paths of every GO-term from itself to the root node.
        file = open('./gos/paths_go_'+anc+'.tsv', 'w')
        #ind: GO-term.
        for ind in child_parent.index:
            #Split the string of parents up into arrays.
            neighbors=child_parent.loc[ind].values[0].split(',')
            #For every parent, their paths to roots are found by depth-first traversal.
            for start in neighbors:
                if start in child_parent.index:
                    traces = real_dfs(child_parent, start)
                    #Of course, every GO-term itself must be included as the starting point. 
                    #At this stage, the duplicated occurrences of GO-term does not matter.
                    #They will be taken cared by the later steps.
                    tmp_str=ind+','
                    for t in traces:
                        tmp_str += t + ','
                    if len(tmp_str) != 0:
                        file.write(ind + '\t'+ tmp_str[:-1]+'\n')
                    else:
                        file.write(ind + '\t'+'\n')


        file.close()

        #Reorganize the TSV-file to the desired format by reading and writing.
        paths=pd.read_csv('./gos/paths_go_'+anc+'.tsv', sep='\t', index_col=0, header=None)
        paths.to_csv('./gos/paths_go_'+anc+'.tsv', sep='\t')
        paths=pd.read_csv('./gos/paths_go_'+anc+'.tsv', sep='\t', index_col=0)


        # concatenate the paths for each GO-terms.
        paths = paths.reset_index().astype(str)
        paths['1']=paths.groupby(['0'])['1'].transform(lambda x : ','.join(x.values))

        # drop duplicate data
        paths = paths.drop_duplicates().set_index('0')

        #Remove the repetitive edges of parent and child go-terms every gene.
        paths['1']=paths['1'].apply(lambda x: ','.join(np.unique(x.split(','))) )
        paths.to_csv('./gos/paths_go_'+anc+'.tsv', sep='\t')


def paths(file_name):
    '''
    A help function that send out every node from paths of each GO-term.
    yield is an efficient alternative of return statement.
    '''
    paths_prob=pd.read_csv(file_name, sep='\t', index_col=0)
    for item in paths_prob['path'].values:
        nodes=item.split(',')
        for n in nodes:
            yield n    
    
def go_counts():
    '''
    For every ontology, the occurrences of each GO-terms are counted.
    '''
    for anc in ['bp', 'cc', 'mf']:
        d = collections.Counter(paths('./gos/paths_prob_'+anc+'.tsv'))
        df_counts = pd.DataFrame.from_dict(d, orient='index')
        df_counts.to_csv('./gos/go_counts_'+anc+'.tsv', sep='\t')



def extract_gos_uniprot():
    '''
    A function that extracts the GO-terms of each gene in the Uniprot database. 
    '''
    evidences = ['EXP', 'IDA','IMP', 'IPI', 'IEP', 'IGI', 'IC', 'TAS']
    #Because of the huge size of the GAF-file, it needs to be imported in chunks.
    for uniprot in pd.read_csv('./gos/goa_uniprot_all.gaf', sep='\t', header=None, skiprows=9, chunksize=50000):
        uniprot.columns=np.arange(len(uniprot.columns)).tolist()
        uniprot=uniprot[[1, 4, 6]]
        uniprot=uniprot.loc[uniprot[6].isin(evidences)][[1, 4]].drop_duplicates().dropna()  
        uniprot.to_csv('./gos/uniprot_extracted.tsv', sep='\t', mode='a', header=False)
            

def paths_prob():
    '''
    This function reorganize the extracted GO-terms of each gene in Uniprot database.
    This will produce a table with every gene and its paths of DFS and also the assigned GO-terms.
    Why does this table need to be produced? 
    Because computing the Yang-Clark and Schlicker similarity needs their paths and assigned GO-terms, respectively.
    '''
    gos_uniprot=pd.read_csv('./gos/uniprot_extracted.tsv', sep='\t', header=None, index_col=0)
    gos_uniprot.columns=['id', 'go_term']
    for anc in ['bp', 'cc', 'mf']:
        paths_go=pd.read_csv('./gos/paths_go_'+anc+'.tsv', sep='\t').dropna()
        #This is to ensure that there is no GO-term that is not included in the original OBO-file.
        tmp=gos_uniprot.loc[gos_uniprot['go_term'].isin(paths_go.set_index('0').index.tolist())]

        #Concatenate the GO-terms for every gene in Uniprot.                                                              
        tmp['go_term']=tmp.groupby('id')['go_term'].transform(lambda x: ','.join(x.values))
        tmp['go_term']=tmp['go_term'].apply(lambda x: ','.join(np.unique(str(x).split(',')))).replace({'nan':float('nan')})
        tmp=tmp.dropna().drop_duplicates()    
        paths_prob=paths_go.merge(gos_uniprot, left_on='0', right_on='go_term')[['id', '1']] 
        paths_prob.columns=['id', 'path']
        paths_prob['path']=paths_prob.drop_duplicates().groupby('id')['path'].transform(lambda x : ','.join(x.values))
        paths_prob['path']=paths_prob['path'].apply(lambda x: ','.join(np.unique(str(x).split(',')))).replace({'nan':float('nan')})
        paths_prob=paths_prob.dropna() 
        paths_prob=tmp.merge(paths_prob, on='id').dropna().drop_duplicates()
        paths_prob.to_csv('./gos/paths_prob_'+anc+'.tsv', sep='\t')
   
        
def bayes_prob_table():
    '''
    Construct a table of the joint counts of all parents and the joint counts of this GO-term with all parents.
    This is one of the three table that will be combined later.
    '''
    for anc in ['bp', 'cc', 'mf']:
        child_parent=pd.read_csv('./gos/child_parent_'+anc+'.tsv', sep='\t', index_col=0).dropna()
        paths_prob=pd.read_csv('./gos/paths_prob_'+anc+'.tsv', sep='\t', index_col=0).dropna()['path'].values.tolist()
        
        #For each child_parent relationship, determine its joint probability with and without child.
        #Joint probability can be computed as the division of two terms:
        #The count with co-occurrences of all parents with and without child.
        collector = []
        for c in child_parent.index:
            cp = child_parent.loc[c].values[0]
            count=0
            count_child=0
            if cp != '':
                count=np.sum(list(map(lambda x: float(all_in(cp.split(','), x.split(','))), paths_prob)))
                count_child=np.sum(list(map(lambda x: float(all_in((c+','+cp).split(','), x.split(','))), paths_prob))) 
     
            collector.append([c, count, count_child])
        bayes=pd.DataFrame(collector, columns=['child_node','joint_count', 'joint_count_child'])
        bayes['cond_prob']=bayes['joint_count_child']/bayes['joint_count']
        bayes.to_csv('./gos/bayes_prob_table_'+anc+'.tsv', sep='\t')


def all_in(query, target):
    '''
    It is a help function of bayes_prob_table().
    It determines if all the queries array(list) are in the target array(list).
    :param query: An array served as query.
    :param target: An array to be determined whether all queries are inside.
    '''
    truth=True
    for q in query:
        if not(q in target):
            truth = False
            break
    
    return truth

def merge3table():
    '''
    Merge bayes_prob_table_xx.tsv, paths_go_xx.tsv and go_counts_xx.tsv together.
    '''
    for anc, go in zip(['bp', 'cc', 'mf'], ['GO:0008150', 'GO:0005575', 'GO:0003674']):
        bayes=pd.read_csv('./gos/bayes_prob_table_'+anc+'.tsv', sep='\t', index_col=0)
        go_counts=pd.read_csv('./gos/go_counts_'+anc+'.tsv', sep='\t')
        paths_go=pd.read_csv('./gos/paths_go_'+anc+'.tsv', sep='\t')    
        go_counts.columns=['node', 'count_node']
        paths_go.columns=['node', 'path']
        bayes=bayes.rename(columns={'child_node': 'node'})
        probs=paths_go.merge(go_counts.merge(bayes, on='node'), on='node')
        probs['count_ancestor']=go_counts.set_index('node').loc[go,'count_node']
        probs['rel_freq']=probs['count_node']/probs['count_ancestor']
        probs.to_csv('./gos/probs_'+anc+'.tsv', sep='\t')



     
def main():
    '''
    Main function calls.
    '''
    #Just to compute the total time.
    t1=time.time()
    child_parent() 
    paths_go()
    extract_gos_uniprot()
    paths_prob()
    bayes_prob_table()
    go_counts()
    merge3table()
    t2=time.time()
    print(t2-t1)
main()

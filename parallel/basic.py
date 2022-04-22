import DSGRN
from DSGRN import *
from multiprocessing import Pool
import itertools
import networkx as nx
import re

import DSGRN
from DSGRN import *
import sys, os

# Disable printing
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def get_network_string(edges, bool):
    """
    edges: a tuple of edges
    bool: a list of 0,1's depicting if edge in first item is repressing (0) or activating (1). 
    
    Example:  edges = (('Hb', 'Gt'),
                ('Hb', 'Kr'),
                ('Hb', 'Kni'),
                ('Gt', 'Hb'),
                ('Gt', 'Kr'),
                ('Gt', 'Kni'),
                ('Kr', 'Hb'),
                ('Kr', 'Gt'))

              bool = [0, 0, 1, 1, 0, 0, 1, 0]
    then the edge from 'Hb' to 'Gt' is repressing while the 'Kr' to 'Hb' edge is activating.

    returns: string for use with DSGRN network input.
    """

    net_dict = {'Hb': [], 'Gt': [], 'Kr': [], 'Kni':[]}

    for i in edges:
        index = edges.index(i)
        net_dict[i[1]] += [(bool[index], i[0])]
        
    new = {}
    for node in net_dict:
        d = collections.defaultdict(list)
        act_str = ''
        rep_str = ''
        for k, *v in net_dict[node]:
            d[k].append(v)
        
        for edge_type in list(d.items()):
            if edge_type[0] == 0:
                rep_str = ''
                for i in edge_type[1]:
                    rep_str += '(~' + i[0] + ')'

            if edge_type[0] == 1:
                act_str = '('
                for i in edge_type[1]:
                    if edge_type[1].index(i) == 0:
                        act_str +=  i[0]
                    else:
                        act_str +=  '+' + i[0]
                act_str += ')'
        new[node] = act_str + rep_str

    return '"""Hb : ' + new['Hb'] + '\n' + 'Gt : ' + new['Gt'] + '\n' + 'Kr : ' + new['Kr'] + '\n' + 'Kni : ' + new['Kni'] + '"""'

def get_FG_layer(graph, current_layer):
    new_layer = []
    while current_layer != []:
        node = current_layer.pop()
        for edge in graph.edges:
            if node == edge[0]:
                new_layer.append(edge[-1])
    return list(set(new_layer))

def get_FG_layer_hex(graph, current_layer):

    FG_layer = [graph.data[node][0] for node in current_layer]
    
    return FG_layer

def get_hex_FG(database, gene):
    
    single_gene_query = SingleGeneQuery(database, gene)
    graph = single_gene_query(2)

    FG = {}
    layer = 0
    current_layer = [0]
    FG[layer] = get_FG_layer_hex(graph, current_layer)
    new_layer = get_FG_layer(graph,current_layer)

    while new_layer !=[]:
        layer +=1
        current_layer = new_layer
        FG[layer] = get_FG_layer_hex(graph, current_layer)
        new_layer = get_FG_layer(graph,current_layer)

    return FG

def get_Hb_Kni_list(database):

    with HiddenPrints():
        Hb_list = {}
        Hb = get_hex_FG(database, 'Hb')
        for key in reversed(list(Hb.keys())):
            Hb_list[len(Hb)-1-key] = Hb[key]  #all code requires Hb_list to be in reverse order, as it simplfies code.

        Kni_list = get_hex_FG(database, 'Kni')

    return Hb_list, Kni_list

def get_number_out_edges_from_string(string, out_edges = {'Hb': 0, 'Gt': 0, 'Kr': 0, 'Kni': 0} ):

    for line in re.findall('\w+\s:.*', string, flags=re.MULTILINE):
        key = re.findall('^\w+', line)[0]
        out_edges[key] = re.findall(key, string, flags=re.MULTILINE).count(key)-1
    
    return out_edges

def convert_FP_list_2_FP_str(lst):
    '''
    lst: list of 4 integerts depicting how many thresholds Hb, Gt, Kr and Kni are above (IN THIS ORDER)
    return: The FP string. I.e., input [0,1,2,1] and return FP { 0, 1, 2, 1} string.
    '''
    return 'FP { ' + str(lst[0]) + ', ' + str(lst[1]) + ', ' + str(lst[2]) + ', ' + str(lst[3]) + ' }'

def get_region_head(Phenotype_pattern, num_out_thresh):

    region_head = {}
    for r in Phenotype_pattern:
        region_head[r] = []
        position = {0:[], 1:[], 2:[], 3:[]}
        
        for i in range(4):
            if Phenotype_pattern[r][i] == '*':
                if num_out_thresh[i]>1:
                    position[i] += [n for n in range(1,num_out_thresh[i])]
                if num_out_thresh[i] == 1:
                    position[i] = [0,1]
            if Phenotype_pattern[r][i] == 'H':
                position[i] += [num_out_thresh[i]]
            if Phenotype_pattern[r][i] == 'L':
                position[i] += [0]

        region_FP = [(a,b,c,d) for a in position[0] for b in position[1] for c in position[2] for d in position[3] ]
        for FP in region_FP:
            #region_head[r] += [convert_FP_list_2_FP_str(FP)]
            region_head[r] += [FP]
    return region_head

def get_FP_Poset(out_edges):
    '''
    out_edges: dictinary with network nodes as keys and number of out edges as values.
    returns: the fixed point pattern between layers we need any path to follow.
    '''
    ### Even if not given in order, this way makes sure we have a particular order
    Phenotype_pattern = {1: '*HLL', 2: 'H*LL', 3: 'HL*L', 4: '*LHL', 5: 'LLH*', 6: 'LL*H', 7: 'L*LH', 8:'LHL*'}
    num_out_thresh = []
    num_out_thresh.append(out_edges['Hb'])
    num_out_thresh.append(out_edges['Gt'])
    num_out_thresh.append(out_edges['Kr'])
    num_out_thresh.append(out_edges['Kni'])

    region_head = get_region_head(Phenotype_pattern, num_out_thresh)
    FP_Regions = {}
    FP_Poset_graph = nx.DiGraph()
    for r in range(1,8):
        FP_Regions[r] = [convert_FP_list_2_FP_str(lst) for lst in region_head[r].copy()]

        temp = nx.DiGraph()
        for n in region_head[r]:
            temp.add_node(n)
        for n in region_head[r+1]:
            temp.add_node(n)

        position = {0:[], 1:[], 2:[], 3:[]}
        inc_dec = {}

        for i in range(4):
            if Phenotype_pattern[r][i] == '*':
                if Phenotype_pattern[r+1][i] == 'H':
                    inc_dec[i] = 'inc'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(1,num_out_thresh[i]+1)]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]
                if Phenotype_pattern[r+1][i] == 'L':
                    inc_dec[i] = 'dec'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(0,num_out_thresh[i])]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]

            if Phenotype_pattern[r+1][i] == '*':
                if Phenotype_pattern[r][i] == 'H':
                    inc_dec[i] = 'dec'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(1,num_out_thresh[i]+1)]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]
                if Phenotype_pattern[r][i] == 'L':
                    inc_dec[i] = 'inc'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(0,num_out_thresh[i])]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]
            if Phenotype_pattern[r][i] == Phenotype_pattern[r+1][i]:
                inc_dec[i] = True
                position[i] += [region_head[r][0][i]]
        
        region_T = [(a,b,c,d) for a in position[0] for b in position[1] for c in position[2] for d in position[3]]  
        
        for t in region_T:
            temp.add_node(t)
        #print(region_T)
        #print(inc_dec)
        FP_Regions[r] += [convert_FP_list_2_FP_str(i) for i in region_T.copy() if i not in region_head[r+1] and i not in region_head[r]]
        for s in temp:
            for t in temp:
                if s != t:
                    edge = True
                    for i in range(4):
                        if inc_dec[i] == True and s[i] != t[i]:
                            edge = False
                        if inc_dec[i] == 'inc' and s[i] > t[i]:
                            edge = False
                        if inc_dec[i] == 'dec' and s[i] < t[i]:
                            edge = False
                    if edge == True:
                        temp.add_edge(s,t)
                        #print(s,t)
        #print(temp.edges())
        for e in temp.edges():
            FP_Poset_graph.add_edge(e[0], e[1])
    mapping = {}
    for node in FP_Poset_graph:
        mapping[node] = convert_FP_list_2_FP_str(node)

    FP_Poset_graph = nx.relabel_nodes(FP_Poset_graph, mapping)

    FP_Poset = nx.to_dict_of_lists(FP_Poset_graph, nodelist=FP_Poset_graph.nodes())
    FP_Regions[8] = [convert_FP_list_2_FP_str(lst) for lst in region_head[8].copy()]

    return FP_Poset, FP_Regions

def single_param(tup):
    params, s = tup
    database_filename, FP_keep, Hb, Kni = params

    database = Database(database_filename)
    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)

    MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(s))
    MGI = MGI_result.fetchone()[0]
    FP_results = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
    if len(FP_results) == 1 and FP_results[0][0:2] == 'FP':
        if set(FP_results).intersection(set(FP_keep)):
            sHb = Hb[((((pg.parameter(s)).logic())[0]).stringify())[6:-2]]
            sKni = Kni[((((pg.parameter(s)).logic())[3]).stringify())[6:-2]]
            for t in list(pg.adjacencies(s, 'codim1')):
                MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(t))
                MGI = MGI_result.fetchone()[0]
                FP_resultt = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
                if len(FP_resultt) == 1 and FP_resultt[0][0:2] == 'FP':
                    if set(FP_resultt).intersection(set(FP_keep)):
                        tHb = Hb[((((pg.parameter(t)).logic())[0]).stringify())[6:-2]]
                        tKni = Kni[((((pg.parameter(t)).logic())[3]).stringify())[6:-2]] 
                        if sHb+1 == tHb and sKni == tKni:
                            print(((sHb, sKni, s), (tHb, tKni, t)), flush=True)
                        elif sHb == tHb and sKni+1 == tKni:
                            print(((sHb, sKni, s), (tHb, tKni, t)), flush=True)
                        elif sHb == tHb and sKni == tKni:
                            print(((sHb, sKni, s), (tHb, tKni, t)), flush=True)

if __name__ == '__main__':
    database = Database("/home/elizabeth/Desktop/HPC/test.db")
    pg = DSGRN.ParameterGraph(database.network)

    with open( "/home/elizabeth/Desktop/HPC/test.txt","r") as f:
        network = f.read()
    out_edges = get_number_out_edges_from_string(network)
    database_filename = "/home/elizabeth/Desktop/HPC/test.db"
    FP_Poset = get_FP_Poset(out_edges)[0]
    FP_keep = [node for node in FP_Poset.keys()]

    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb = {}
    for i in Hb_list:
        for j in Hb_list[i]:
            Hb[j] = i
    Kni = {}
    for i in Kni_list:
        for j in Kni_list[i]:
            Kni[j] = i

    with Pool(4) as p:
        p.map(single_param, zip(itertools.repeat((database_filename, FP_keep, Hb, Kni)), range(0,pg.size())))
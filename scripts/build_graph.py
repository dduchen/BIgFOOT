#!/usr/bin/python3

import sys
import os.path
from graph_tool.all import *
import time
import math
import subprocess
import argparse
import re


__author__ = "Jasmijn Baaijens"
__license__ = "MIT"


usage = """
Build a variation graph from a set of contigs using a partial order alignment.
"""


def main():
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('contigfile', type=str, help='The file containing the contigs, format is .fasta')
    parser.add_argument('-p', type=str, default='poa_graph.txt', help='Store POA output graph here')
    parser.add_argument('-s', type=str, help='The alignment scoring matrix; by default the blosum80 matrix is used.')
    parser.add_argument('-g', type=str, default='contig_graph.gt', help='Store the variation graph here')
    parser.add_argument('-o', type=str, default='graphviz.png', help='PNG file to store the visualisation of the graph')
    parser.add_argument('--skip_poa', action='store_true', help = 'Skip re-doing POA alignment; assumes there is an existing alignment file')
    parser.add_argument('--visualize', action='store_true', help='Option to visualize the graph')

    args = parser.parse_args()

    # Perform a POA and make the initial graph
    if args.s:
        blosumfile = args.s
    else:
        blosumfile = os.path.dirname(os.path.abspath(__file__)) + "/../blosum80.mat"
    g = make_graph(args.contigfile, args.p, blosumfile, args.skip_poa)

    # Check if graph is cyclic
    is_cyclic = cyclic(g)
    if is_cyclic:
        print("WARNING: initial graph contains cycles!!!")
    else:
        print("The initial graph is acyclic :)")

    # Compress the inital graph
    g = compress_graph(g)
    #write_gfa(g, "intermediate_graph.gfa")

    # # Remove duplicate nodes
    # removed_count = 1
    # while removed_count > 0:
    #     # Deduplicate
    #     g, removed_count = merge_duplicate_nodes(g)
    #     print("removed {} duplicate nodes".format(removed_count))
    #     # Compress
    #     g = compress_graph(g)
    #     # Check if graph is cyclic
    #     is_cyclic = cyclic(g)
    #     if is_cyclic:
    #         print("WARNING: graph contains cycles!!!")
    #     else:
    #         print("The graph is acyclic :)")

    # Check if graph is cyclic
    is_cyclic = cyclic(g)
    if is_cyclic:
        print("WARNING: final graph contains cycles!!!")
    else:
        print("The final graph is acyclic :)")

    # Save the graph to the graphfile & GFA
    g.save(args.g)
    gfa_file = args.g.rstrip('.gt') + '.gfa'
    write_gfa(g, gfa_file, write_paths=True)

    # Graph visualization
    gprops = {'K': 0.5, 'ordering': 'out', 'rank': 'source', 'rankdir': 'LR', 'ratio': 'auto'}
    vprops_seq = {'label': g.vp.seq, 'xlabel':g.vp.contigs,'shape': 'rect', 'fillcolor': g.vp.col,'margin':0}
    eprops = {'arrowsize': 2.0, 'dir': 'forward'}
    if args.visualize:
        graphviz_draw(g, output= args.o, layout = 'dot', size = (500,500), vsize = 1, overlap=False, gprops = gprops, vprops= vprops_seq, eprops = eprops)

    print('Runtime is', time.clock())

    return


def write_gfa(g, gfa_file, write_paths=False):
    """
    Write graph-tool graph to a gfa format storing only nodes and edges.
    Returns nothing.
    """
    contigs = {}
    with open(gfa_file, 'w') as file:    #Write the nodes and edges in GFA format
        file.write("H\tVN:Z:1.0")
        for v in g.vertices():
            file.write("\nS\t%d\t%s" % (int(v)+1, g.vp.seq[v]))
            for w in v.out_neighbors():
                file.write("\nL\t%d\t+\t%d\t+\t0M" % (int(v)+1, int(w)+1))
            for k in g.vp.contigs[v]:
                k = int(k)
                if k in contigs.keys():
                    contigs[k].append(v)
                else:
                    contigs[k] = []
                    contigs[k].append(v)

        if write_paths:
            for k in range(len(contigs.keys())):
                path = ''
                match = ''
                for v in contigs[k]:
                    path += str(int(v)+1)+'+,'
                    match += str(len(g.vp.seq[v]))+'M,'
                path = path.rstrip(',')
                match = match.rstrip(',')
                file.write("\nP\tp%d\t%s\t%s" % (k, path, match))
                # file.write("\nP\tp%d\t%s\t*" % (k, path))
    return



def make_graph(contigfile, poafile, scoring_matrix, skip_poa):
    """
    Run POA on contig file and POA scoring matrix, then build initial contig
    variation graph from POA output. Returns graph in graph-tool format.
    """
    #Perform a progressive multiple sequence alignment
    if not skip_poa:
        subprocess.check_call(['poa', '-read_fasta', contigfile, '-po', poafile,
                    scoring_matrix, '-do_progressive',
                    '-preserve_seqorder'])

    g = Graph(directed=True) #Define a graph with its vertex properties
    vprop = g.new_vertex_property('string')
    g.vp.seq = vprop
    vprop = g.new_vertex_property('vector<string>')
    g.vp.contigs = vprop

    f = open(poafile, 'r')  #Open the file with the POA results

    c_start = []
    v_start = []
    header = True
    for line in f:
        line = line.strip().split(':')
        base = line[0]
        if base in ['A','a','T','t','G','g','C','c','N']:
            # POA format: <node_id>:<node_info> where node info is identified by
            # residue labels 'L' = incoming edges
            #                'S' = sequences stored in node
            #                'A' = index of next node in same align ring
            header = False
            rest = line[1]
            # store node in graph
            v = g.add_vertex()
            g.vp.seq[v] = base.upper()
            # process node info: add edges and store contigs
            info = list(filter(None, re.split(r'(\d+)', rest)))
            assert len(info)%2 == 0
            i = 0
            while i < len(info)-1:
                if info[i] == 'L':
                    # store edge
                    node = info[i+1]
                    g.add_edge(node, v)
                elif info[i] == 'S':
                    # store contig
                    contig = info[i+1]
                    g.vp.contigs[v].append(contig)
                else:
                    assert info[i] == 'A'
                i += 2
        elif not header:
            # skip header
            print("ERROR: base not recognized. Exiting.")
            sys.exit(1)
    f.close()
    return g


def merge_duplicate_nodes(g):
    """
    Remove duplicate nodes (i.e. parallel nodes with identical sequence) from
    the input graph. WARNING: POA adds such nodes to make the graph acyclic,
    hence merging causes cycles.
    Returns the updated graph (gt-format) and the number of removed nodes.
    """
    del_list = []
    for v in list(g.vertices()):
        if v.out_degree() > 1:
            base_to_node = {}
            for w in sorted(list(v.out_neighbors())):
                if w in del_list:
                    continue
                b = g.vp.seq[w]
                if b in base_to_node:
                    existing_node = base_to_node[b]
                    contigs1 = list(g.vp.contigs[existing_node])
                    contigs2 = list(g.vp.contigs[w])
                    stop = False
                    for neighbor in list(w.in_neighbors()):
                        if neighbor == existing_node:
                            stop = True
                        if neighbor not in list(existing_node.in_neighbors()):
                            g.add_edge(neighbor, existing_node)
                    for neighbor in list(w.out_neighbors()):
                        if neighbor == existing_node:
                            stop = True
                        if neighbor not in list(existing_node.out_neighbors()):
                            g.add_edge(existing_node, neighbor)
                    if not stop:
                        g.vp.contigs[existing_node] = list(set(contigs1 + contigs2))
                        del_list.append(w)
                else:
                    base_to_node[b] = w

        if v.in_degree() > 1:
            base_to_node = {}
            for w in sorted(list(v.in_neighbors())):
                if w in del_list:
                    continue
                b = g.vp.seq[w]
                if b in base_to_node:
                    existing_node = base_to_node[b]
                    contigs1 = list(g.vp.contigs[existing_node])
                    contigs2 = list(g.vp.contigs[w])
                    stop = False
                    for neighbor in list(w.in_neighbors()):
                        if neighbor == existing_node:
                            stop = True
                        if neighbor not in list(existing_node.in_neighbors()):
                            g.add_edge(neighbor, existing_node)
                    for neighbor in list(w.out_neighbors()):
                        if neighbor == existing_node:
                            stop = True
                        if neighbor not in list(existing_node.out_neighbors()):
                            g.add_edge(existing_node, neighbor)
                    if not stop:
                        g.vp.contigs[existing_node] = list(set(contigs1 + contigs2))
                        del_list.append(w)
                else:
                    base_to_node[b] = w

    # now remove all selected vertices
    g.remove_vertex(del_list)
    return g, len(del_list)


def compress_graph(g):
    """
    Compress the variation graph: merge any non-branching paths into a single
    node and update all node properties accordingly. Returns the updated graph.
    """
    #(g, v_start) = graph
    vprop = g.new_vertex_property('string') #Make a new vertex property: color
    g.vp.col = vprop
    del_list = []
    for v in list(g.vertices()):
        g.vp.col[v] = 'white'                   #Colour the vertex white

        if v.out_degree() == 0:                 #If a vertex is not connected
            if v.in_degree() == 0:              #to the graph, then remove it
                # for i in range(len(v_start)):
                #     if v_start[i] >= int(v):
                #         v_start[i] += -1
                del_list.append(v)
        elif v.out_degree() == 1:
            w = list(v.out_neighbors())[0]
            if w.in_degree() == 1:
                contigs_v = list(g.vp.contigs[v])
                contigs_w = list(g.vp.contigs[w])
                g.vp.seq[w] = ''.join([g.vp.seq[v], g.vp.seq[w]])
                g.vp.contigs[w] = list(set(contigs_v + contigs_w))
                g.vp.col[w] = g.vp.col[v]
                v_in_neighbors = list(v.in_neighbors())
                for neighbor in v_in_neighbors:
                    #if neighbor not in list(w.in_neighbors()):
                    g.add_edge(neighbor, w)
                del_list.append(v)

    # now remove all selected vertices
    g.remove_vertex(del_list)

    length_g = 0
    for v in g.vertices():
        length_g += len(g.vp.seq[v])  #Sum the number of letters used
        if v.out_degree() == 0 or v.in_degree() == 0: #If the graph start or
            g.vp.col[v] = 'grey'                      #ends in v, colour it grey
        # if v in v_start:                              #If a contig start in v
        #     g.vp.col[v] = 'grey'                      #also colour it grey

    print('The total sequence length of the graph is',length_g)
    return g



def cyclic(graph):
    """Return True if the directed input graph has a cycle."""
    visited = set()
    path = [object()]
    path_set = set(path)
    stack = [graph.vertices()]
    while stack:
        for v in stack[-1]:
            if v in path_set:
                print(list(path_set))
                print(v)
                return True
            elif v not in visited:
                visited.add(v)
                path.append(v)
                path_set.add(v)
                stack.append(iter(v.out_neighbors()))
                break
        else:
            path_set.remove(path.pop())
            stack.pop()
    return False


if __name__ == '__main__':
    sys.exit(main())

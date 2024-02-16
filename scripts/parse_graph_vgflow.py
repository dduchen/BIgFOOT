#!/usr/bin/env python3

import sys, os
import json
#sys.path.append('/Users/dylanduchen/tools/vg-flow/scripts/')
sys.path.append('/home/dd392/tools/vg-flow/scripts/')
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort
import subprocess
import argparse
# from tqdm import tqdm # progress tracker
import time
from graph_functions import has_simple_paths, reorder_graph, compress_graph#, get_node_abundances
from graph_functions import write_gfa, read_gfa, cyclic, check_parallel_nodes
from graph_functions import filter_edges


__author__ = "Adaptations by Dylan Duchen - Original code by Jasmijn Baaijens"
__license__ = "MIT"

usage = "Parse contig variation graph and compute node abundances."


def main():
    parser = argparse.ArgumentParser(prog='parse_graph_vgflow.py', description=usage)
    parser.add_argument("--sample", default="SRR14814127", help="sample ID and core component of filenames")
    parser.add_argument('-m', '--min_edge_ab', dest='min_edge_ab', type=int, default=0, help="Minimal number of reads spanning an edge; edges with abundance below this threshold are removed from the contig variation graph")
    parser.add_argument('--filter_branches_only', action='store_true', default=True, help="Only apply edge abundance filter to branching edges")
    args = parser.parse_args()


    filepath = os.path.dirname(os.path.abspath(__file__))
    graph_name = args.sample
    gfa_file = "{}.gfa".format(graph_name)
    aln_name = "{}.aln".format(graph_name)
    aln_file = "{}.json".format(aln_name)

    contig_graph, paths, node_dict = read_gfa(gfa_file)
    print("Computing node/edge abundances...",)

    node_abundances, vertex_pair_abundances = get_node_abundances(
        contig_graph, node_dict, aln_file)

    # write edge abundances to file
    edge_abundances = {}
    with open('{}.edge_abundances.txt'.format(graph_name), 'w') as f:
        for e in contig_graph.edges():
            v1 = e.source()
            v2 = e.target()
            try:
                edge_ab = vertex_pair_abundances[(v1, v2)]
            except KeyError as err:
                edge_ab = 0
            edge_abundances[e] = edge_ab
            f.write("{},{}\t{}\n".format(v1, v2, edge_ab))

    # remove edges of low abundance and break contigs accordingly
    write_gfa(contig_graph, '{}.tmp1.gfa'.format(graph_name), paths=paths)
    contig_graph, paths = filter_edges(contig_graph, paths, edge_abundances,
                                       args.min_edge_ab,
                                       args.filter_branches_only)

    # merge simple paths and compute average abundance counts
    print("Compressing graph...",)
    write_gfa(contig_graph, '{}.tmp2.gfa'.format(graph_name), paths=paths)
    final_graph, final_paths = compress_graph(
            contig_graph, paths, node_abundances)
    V = len(list(final_graph.vertices()))
    E = len(list(final_graph.edges()))
    print("Final graph size: {} nodes and {} edges".format(V, E))
    # write_gfa(final_graph, "compressed_graph.gfa", paths=final_paths)

    # reorder nodes in topological order
    if is_DAG(final_graph):
        print("Final graph is acyclic :)")
        top_ordering = topological_sort(final_graph)
        final_graph, final_paths = reorder_graph(
                final_graph, final_paths, top_ordering)

    # write final graph and node abundances to file
    write_gfa(final_graph, '{}.final.gfa'.format(graph_name), paths=final_paths)
    final_graph.save('{}.final.gt'.format(graph_name))
    abundance_file = open('{}.node_abundance.txt'.format(graph_name), 'w')
    for node in final_graph.vertices():
        v = int(node)
        if final_graph.vp.seq[node] == "N":
            # make sure ambiguous bases get low abundance to allow correction
            ab = 0 # also gets ignored during objective evaluation
        else:
            ab = final_graph.vp.ab[node]
        abundance_file.write("{}:{}\n".format(v, ab))

    abundance_file.close()

    print("\ngraph processing completed")
    print()
    return

#####################################################

def get_node_abundances(graph, node_dict, aln_file):
    """
    Reads vg alignments and computes node abundance values and edge abundances
    (connection counts per edge).
    Returns node abundance list and connection count dict.
    """
    node_dict_new2old = {v_new : v_old for v_old, v_new in node_dict.items()}
    nodes = {}
    for node in graph.vertices():
        # update node ID, because alignments are to GFA graph which has 1-based
        # node IDs
        node_id = node_dict_new2old[int(node)]
        seq = graph.vp.seq[node]
        nodes[node_id] = seq
    bases_per_node = {} # map node IDs to read alignments
    for node in nodes:
        bases_per_node[node] = 0
    # count node connections used by read alignments
    connection_counts = {}
    print("Processing alignments...")
    with open(aln_file, 'r') as aln_json:
        # for line in tqdm(aln_json):
        for line in aln_json:
            aln = json.loads(line)
            seq_id = aln["name"]
            seq = aln["sequence"]
            try:
                path = aln["path"]
                mapping = path["mapping"]
            except KeyError:
                # read unmapped
                continue
            offset = 0
            seq_aln_len = 0
            node0 = ""
            node1 = ""
            node2 = ""
            for node_info in mapping:
                position = node_info["position"]
                node_id = position["node_id"]
                node2 = node_dict[node_id]
                node_seq = nodes[node_id]
                node_len = len(node_seq)
                aln_len = 0
                try:
                    offset = int(position["offset"])
                except KeyError:
                    offset = 0
                try:
                    is_reverse = bool(position["is_reverse"])
                except KeyError:
                    is_reverse = False
                edit = node_info["edit"]
                for aln_piece in edit:
                    try:
                        from_len = int(aln_piece["from_length"])
                    except KeyError:
                        from_len = 0
                    try:
                        to_len = int(aln_piece["to_length"])
                    except KeyError:
                        to_len = 0
                    aln_len += min(from_len, to_len)
                # bases_per_node[node_id] += node_len - offset
                bases_per_node[node_id] += aln_len
                #node_to_seq[node_id].append(seq_id)
                seq_aln_len += aln_len
                if node1 != "":
                    # store direct connection
                    if not is_reverse:
                        # assert node2 in graph.vertex(node1).out_neighbors()
                        try:
                            connection_counts[(node1, node2)] += 1
                        except KeyError:
                            connection_counts[(node1, node2)] = 1
                    else:
                        # assert node1 in graph.vertex(node2).out_neighbors()
                        try:
                            connection_counts[(node2, node1)] += 1
                        except KeyError:
                            connection_counts[(node2, node1)] = 1
                node0 = node1
                node1 = node2
    # TODO: check for double alignments and weight them accordingly
    # for node, seq_list in node_to_seq.items():
    #     for seq in seq_list:
    #         if seq_list.count(seq) > 2:
    #     #if len(set(seq_list)) != len(seq_list):
    #             print("DUPLICATE READ {} IN SEQ LIST OF NODE {}".format(seq, node))
    #             print(seq_list)
    #             sys.exit(1)
    # node_abundance_file = "node_abundance.txt"
    node_abundance_list = {}
    print("Computing node abundance rates...")
    # with open(node_abundance_file, 'w') as f:
    for node, seq in nodes.items():
        node_len = len(seq)
        aligned_len = bases_per_node[node]
        if node_len > 0:
            node_abundance = aligned_len / node_len
        else:
            print("Node length 0 for node {} ?!".format(node))
            node_abundance = 0
        node_abundance_list[node_dict[node]] = node_abundance
        # print(node, ":", node_abundance)
        # f.write("{0}:{1}\n".format(node-1, node_abundance))
    return node_abundance_list, connection_counts


if __name__ == '__main__':
    sys.exit(main())
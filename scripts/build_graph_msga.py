#!/usr/bin/env python3

import sys, os
import json
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort
import subprocess
import argparse
# from tqdm import tqdm # progress tracker
import time

from graph_functions import has_simple_paths, reorder_graph, compress_graph
from graph_functions import write_gfa, read_gfa, cyclic, check_parallel_nodes
from graph_functions import filter_edges

__author__ = "Jasmijn Baaijens"
__license__ = "MIT"

usage = "Build contig variation graph using vg msga and compute node abundances."


def main():
    parser = argparse.ArgumentParser(prog='build_graph_msga.py', description=usage)
    parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads in fastq format')
    parser.add_argument('-r', '--reverse', dest='reverse', type=str, help='Reverse reads in fastq format')
    parser.add_argument('-c', '--contigs', dest='contigs', type=str, required=True, help='Input contigs in fastq format')
    parser.add_argument('-vg', '--vg_path', dest='vg', type=str, required=True, help="Path to vg executable")
    parser.add_argument('-t', '--threads', default=8, help="Set number of threads used for vg.")
    parser.add_argument('-w', '--msga_w', dest='msga_w', default=128, help="Set alignment band width parameter for vg msga.")
    parser.add_argument('-m', '--min_edge_ab', dest='min_edge_ab', type=int, default=0, help="Minimal number of reads spanning an edge; edges with abundance below this threshold are removed from the contig variation graph")
    parser.add_argument('--filter_branches_only', action='store_true', help="Only apply edge abundance filter to branching edges")
    parser.add_argument('--reuse_graph', dest='reuse_graph', action='store_true', help="Use contig variation graph stored in 'contig_graph.norm.vg'")
    parser.add_argument('--reuse_index', dest='reuse_index', action='store_true', help="Use contig variation graph indexes stored in 'contig_graph.norm.xg' and 'contig_graph.norm.gcsa'")
    parser.add_argument('--reuse_aln', dest='reuse_aln', action='store_true', help="Use read alignments in 'contig_graph.norm.aln.filtered.json'")
    parser.add_argument('--quick', action='store_true', help="sort contigs in quick mode; use this option if sorting contigs is a bottleneck")
    args = parser.parse_args()

    global_t1_start = time.perf_counter()
    global_t2_start = time.process_time()

    filepath = os.path.dirname(os.path.abspath(__file__))
    graph_name = "contig_graph"
    gfa_file = "{}.norm.gfa".format(graph_name)
    aln_name = "{}.norm.aln".format(graph_name)
    aln_file = "{}.filtered.json".format(aln_name)
    vg = args.vg
    
    if not args.reuse_graph:
        t1_start = time.perf_counter()
        t2_start = time.process_time()
        # sort contigs
        sorted_contigs = "sorted_contigs.fasta"
        if args.quick:
            subprocess.check_call(
                "python {}/sort_contigs.py --quick -t {} {} {}".format(
                    filepath, args.threads, args.contigs, sorted_contigs),
                shell=True
            )
        else:
            subprocess.check_call(
                "python {}/sort_contigs.py -t {} {} {}".format(filepath,
                                                               args.threads,
                                                               args.contigs,
                                                               sorted_contigs),
                shell=True
            )
        with open(sorted_contigs, 'r') as f:
            line = f.readline()
            base_seq = line.lstrip('>').rstrip()

        # build initial graph
        print("Building initial graph")
        subprocess.check_call("mkdir -p tmp", shell=True)
        subprocess.check_call("export TMPDIR=tmp && " + vg +
            " msga -f {} -b {} -t {} -a -w {} > {}.vg".format(
                sorted_contigs, base_seq, args.threads, args.msga_w, graph_name),
            shell=True
        )

        # normalize graph
        print("Normalizing graph")
        subprocess.check_call( vg +
            " mod -n {0}.vg | {1} mod -X 32 - | {1} ids -c - > {0}.norm.vg".format(
                graph_name, vg),
            shell=True
        )

        # convert vg to gfa
        print("Convert graph to GFA")
        subprocess.check_call( vg +
            " view {}.norm.vg > {}".format(graph_name, gfa_file),
            shell=True
        )

        t1_stop = time.perf_counter()
        t2_stop = time.process_time()
        print("Variation graph is ready")
        print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
        print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
        print()

    if not args.reuse_index:
        # index graph
        t1_start = time.perf_counter()
        t2_start = time.process_time()
        print("Building index for graph")
        subprocess.check_call("mkdir -p tmp", shell=True)
        subprocess.check_call( "export TMPDIR=tmp && " + vg +
            " index -x {0}.norm.xg -g {0}.norm.gcsa -k 16 -X 2 -Z 500 -t {1} {0}.norm.vg".format(
                graph_name, args.threads),
            shell=True
        )
        t1_stop = time.perf_counter()
        t2_stop = time.process_time()
        print("Graph index is ready")
        print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
        print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
        print()

    if not args.reuse_aln:
        # align reads to graph
        t1_start = time.perf_counter()
        t2_start = time.process_time()
        print("Mapping reads to graph")
        if args.reverse:
            subprocess.check_call( vg +
                " map -f {0} -f {1} -x {2}.norm.xg -g {2}.norm.gcsa -k 16 -t {4} > {3}.gam".format(
                    args.forward, args.reverse, graph_name, aln_name, args.threads),
                shell=True
            )
        else:
            subprocess.check_call( vg +
                " map -f {0} -x {2}.norm.xg -g {2}.norm.gcsa -k 16 -t {4} > {3}.gam".format(
                    args.forward, args.reverse, graph_name, aln_name, args.threads),
                shell=True
            )

        # filter for primary alignments
        print("Filtering alignments...")
        try:
            subprocess.check_call( vg +
                " filter -r 0.90 -s 2 -fu -t {1} {0}.gam > {0}.filtered.gam".format(
                aln_name, args.threads),
                shell=True)
        except CalledProcessError as e:
            print("WARNING: vg filter failed, continuing with unfiltered reads")
            subprocess.check_call("cp {0}.gam {0}.filtered.gam".format(aln_name))
        t1_stop = time.perf_counter()
        t2_stop = time.process_time()
        print("Alignments ready")
        print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
        print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
        print()

        # output json
        print("Converting gam to json format")
        subprocess.check_call( vg +
            " view -a {0}.filtered.gam > {1}".format(aln_name, aln_file),
            shell=True)

        # clean up vg tmpdir
        subprocess.check_call("rm -rf tmp", shell=True)

    # read graph and check for cycles
    contig_graph, paths, node_dict = read_gfa(gfa_file)
    # # check for parallel sequences in graph
    # check_parallel_nodes(contig_graph)

    # compute node abundances corresponding to the GFA file
    print("Computing node/edge abundances...",)
    node_abundances, vertex_pair_abundances = get_node_abundances(
            contig_graph, node_dict, aln_file)
    # write edge abundances to file
    edge_abundances = {}
    with open('edge_abundances.txt', 'w') as f:
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
    abundance_file = open('node_abundance.txt', 'w')
    for node in final_graph.vertices():
        v = int(node)
        if final_graph.vp.seq[node] == "N":
            # make sure ambiguous bases get low abundance to allow correction
            ab = 0 # also gets ignored during objective evaluation
        else:
            ab = final_graph.vp.ab[node]
        abundance_file.write("{}:{}\n".format(v, ab))
    abundance_file.close()
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    print("\nbuild_graph_msga completed")
    print("Elapsed time: {:.1f} seconds".format(t1_stop-global_t1_start))
    print("CPU process time: {:.1f} seconds".format(t2_stop-global_t2_start))
    print()
    return

######################

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

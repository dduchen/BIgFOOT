#!/usr/bin/python3

import sys
from numpy import array, zeros, multiply, dot, ceil, where, mod
from gurobipy import *
from graph_tool.all import Graph
from graph_tool.topology import is_DAG, all_circuits
import argparse
from tqdm import tqdm # progress tracker
import datetime
import time
import math
from collections import OrderedDict


__author__ = "Jasmijn Baaijens"
__license__ = "MIT"

usage = 'Build strains from contigs in variation graph and estimate their abundances.'


def main():
    parser = argparse.ArgumentParser(prog='python opt_edge.py', description=usage)
    parser.add_argument('abundancefile', type=str, help='Node abundance file')
    parser.add_argument('graphfile', type=str, help='GFA file containing the contig variation graph')
    parser.add_argument('-m', '--min_abundance', dest='m', type=int, required=True, help="Minimal node abundance; nodes below this threshold are considered erroneous and allowed to match any other node.")
    parser.add_argument('-c', '--min_cov', dest='min_cov', type=float, required=True, help='Minimum coverage required per strain')
    parser.add_argument('-d', '--min_depth', dest='min_depth', type=int, default=50, help='Output a list of nodes with sequence depth less than <min_depth>.')
    parser.add_argument('-o', '--out_fasta', dest='fasta', type=str, default='haps.final.fasta', help='Output fasta file with final haplotypes')
    parser.add_argument('-p', '--paths', type=str, help='Read candidate paths from file instead of enumerating them.')
    parser.add_argument('-r', '--reduce_obj', dest='reduce_obj', type=int, default=0, help='Use a reduced objective function by specifying how many times a given combination of paths will be evaluated (reduces runtime and memory usage).')
    parser.add_argument('-t', '--threads', type=int, default=1, help="Set number of threads used for Gurobi.")
    #parser.add_argument('-cb', '--check_branches', action='store_true', help='disable paths that cross branches without contig evidence when available.')
    parser.add_argument('--trim', dest='trim', type=int, default=10, help='number of bases trimmed on either end of contig')
    parser.add_argument('--max_strains', dest='max_strains', type=int, default=0, help='set an upper bound on the number of strains predicted')
    #parser.add_argument('-sp', '--store_paths', action='store_true', help='store all paths and strains in text files.')
    args = parser.parse_args()

    global_t1_start = time.perf_counter()
    global_t2_start = time.process_time()

    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    print()

    # read node abundances from file
    abundance_list = []
    low_ab_nodes = []
    min_ab = args.m
    old_id = -1
    with open(args.abundancefile, 'r') as f:
        for line in f:
            [node_id, abundance] = line.rstrip().split(':')
            # make sure that nodes are ordered incrementally
            new_id = int(node_id)
            if old_id != -1:
                assert new_id == old_id + 1
            old_id = new_id
            # add node abundance to list
            abundance_list.append(float(abundance))
            if float(abundance) < min_ab:
                # store low-abundance nodes
                # convert .gfa node ID (1-based) to .gt node ID (0-based)
                low_ab_nodes.append(int(node_id)-1)

    # read graph from file
    #graph = load_graph(args.graphfile)
    graph, paths = read_gfa(args.graphfile)
    nvert = len(list(graph.vertices()))
    node_lengths = [len(graph.vp.seq[v]) for v in graph.vertices()]
    # count contigs
    contig_IDs = set()
    for contigs in graph.vp.contigs:
        contig_IDs = contig_IDs.union(set(contigs))
    ncontigs = len(contig_IDs)
    print("********* input graph *********")
    print("#vertices = {}".format(nvert))
    print("#edges = {}".format(len(list(graph.edges()))))
    print("#contigs = {}".format(ncontigs))
    print("*******************************\n")

    if args.paths:
        # read paths from file
        vg_paths = []
        with open(args.paths, 'r') as f:
            for line in f:
                if line[0] == '>':
                    continue
                path = [int(x) for x in line.rstrip().split(' ')]
                vg_paths.append(path)
        # construct corresponding sequences
        haps = []
        idx = 0
        with open('haps.fasta', 'w') as f:
            for path in vg_paths:
                seq = ""
                for v in path:
                    seq += graph.vp.seq[graph.vertex(v)]
                haps.append(seq)
                f.write(">{}\n{}\n".format(idx, seq))
                idx += 1
    else:
        t1_start = time.perf_counter()
        t2_start = time.process_time()
        # build all paths consisting of non-conflicting contigs through the graph
        graph, paths, adj_out, start_info = build_adj_info(graph, paths, nvert,
                    low_ab_nodes, args.trim)
        # write contigs and their subpaths to files for error investigation
        write_intermediate_files(graph, paths, adj_out)
        # check if contig overlap graph is cyclic
        cycle_list = cyclic(adj_out)
        if len(cycle_list) > 0:
            # print("contig adjacency graph is cyclic, exiting")
            # sys.exit(1)
            print("contig adjacency graph contains {} cycles.".format(
                    len(cycle_list)))
            print("breaking cycles")
            for cycle in cycle_list:
                print("cycle: ", cycle)
                broken = False
                for v in cycle:
                    if len(adj_out[v]) > 1:
                        print("removing edges {} -> {}".format(v, adj_out[v]))
                        adj_out[v] = []
                        broken = True
                        break
                if not broken:
                    v = cycle[-1]
                    print("removing edges {} -> {}".format(v, adj_out[v]))
                    adj_out[v] = []
        else:
            print("contig adjacency graph is acyclic")

        # check branches: if any contigs bridge the branch, ignore all contigs that
        # start/stop at the internal branch node
        discarded = check_branches2(graph, adj_out)
        # discarded = check_branches(graph, ncontigs, nvert)
        cg_paths = enumerate_contig_paths(graph, adj_out, discarded)
        # print(cg_paths)
        vg_paths, haps = contig_paths_to_nodes(cg_paths, graph, paths, start_info,
                            low_ab_nodes, abundance_list, min_ab)
        # seq_lengths = [len(hap) for hap in haps]
        t1_stop = time.perf_counter()
        t2_stop = time.process_time()
        print("\nPath generation completed")
        print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
        print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))

    skip_nodes = get_extremity_nodes(graph, args.min_depth)
    print("\nFound {} extremity nodes (depth < {})\n".format(len(skip_nodes),
                                                                args.min_depth))

    # if all nodes of a path are skipped we have no constraints -> remove path
    for idx, path in enumerate(vg_paths):
        remove_path = True
        for node in path:
            if node not in skip_nodes:
                remove_path = False
                break
        if remove_path:
            print("Removing path because all nodes are skipped during objective evaluation")
            print("Path:", path)
            paths[idx] = []
    while vg_paths.count([]) > 0:
        vg_paths.remove([])

    # Now solve the minimization problem
    minimization_min_cov = 0
    x, objVal, diff_counter = optimize(abundance_list, nvert, vg_paths,
            skip_nodes, args.reduce_obj, args.max_strains, minimization_min_cov,
            args.min_cov, args.threads)
    # x, objVal, diff_counter = optimize2(abundance_list, nvert, vg_paths,
    #         skip_nodes, args.max_strains, node_lengths, minimization_min_cov,
    #         args.min_cov, args.threads)

    # Analyze optimization output
    final_abundances = process_output(x, diff_counter, haps, args.min_cov, args.fasta)
    build_genome_graph(graph, vg_paths, final_abundances)

    x_stop = time.perf_counter()
    t2_stop = time.process_time()
    print("\noptimize_strains completed")
    print("Elapsed time: {:.1f} seconds".format(t1_stop-global_t1_start))
    print("CPU process time: {:.1f} seconds".format(t2_stop-global_t2_start))
    print()

    return

################################################################################



def read_gfa(gfa_file):
    """
    Reads a graph from a GFA-file and returns graph in gt-format.
    """
    # Define a graph with its vertex properties
    g = Graph(directed=True)
    vprop = g.new_vertex_property('string')
    g.vp.seq = vprop
    vprop = g.new_vertex_property('vector<string>')
    g.vp.contigs = vprop

    # read gfa and add vertices to graph
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'S':
                # vertex
                node_id = int(line[1]) - 1
                v = g.add_vertex()
                seq = line[2].upper()
                if len(seq) == 0:
                    print(line)
                g.vp.seq[v] = seq

    # parse through gfa again to add edges and contig paths to graph
    path_count = 0
    paths = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'L':
                # edge
                assert line[2] == line[4] == "+"
                v1 = int(line[1]) - 1
                v2 = int(line[3]) - 1
                g.add_edge(v1, v2)
            elif line[0] == 'P':
                path_count += 1
                contig_id = line[1]
                path = line[2].split(',')
                unoriented_path = []
                old_ori = ''
                if len(path) == 0:
                    print("ERROR: path of length 0. Exiting.")
                    sys.exit(1)
                for node_info in path:
                    # take care of node orientations
                    if len(node_info) > 0:
                        ori = node_info[-1]
                    else:
                        print("ERROR: no node info in path. Exiting.")
                        sys.exit(1)
                    if old_ori != '':
                        # check if orientations change within the path,
                        # indicating the presence of inversions
                        assert ori == old_ori
                    old_ori = ori
                    node = int(node_info.rstrip(ori)) - 1
                    unoriented_path.append(node)
                    # store contig in node
                    g.vp.contigs[node].append(contig_id)
                assert len(unoriented_path) > 0
                if ori == "+":
                    paths[contig_id] = unoriented_path
                else:
                    paths[contig_id] = unoriented_path[::-1]
    # # save graph in gt format
    # g.save(graph_name + '.gt')
    ordered_paths = OrderedDict(sorted(paths.items(), key=lambda x:x[0]))
    return g, ordered_paths


def build_genome_graph(graph, paths, final_ab):
    """
    Remove nodes that are not traversed by any path of nonzero abundance from
    the variation graph and store haplotype paths, giving the genome variation
    graph. Saves this graph to a gt-file and GFA-file. Returns nothing.
    """
    vprop = graph.new_vertex_property('vector<string>')
    graph.vp.haps = vprop
    del_list = [1 for v in graph.vertices()]
    final_paths = []
    path_num = 0
    for i, path in enumerate(paths):
        if final_ab[i] > 0:
            for v in path:
                del_list[int(v)] = 0
                graph.vp.haps[v].append(str(path_num))
            path_num += 1
            final_paths.append(path)
    # remove non-traversed nodes from graph
    old_to_new_nodes = {}
    i_new = len(del_list) - sum(del_list) - 1
    for i in reversed(range(len(del_list))):
        if del_list[i] == 1:
            graph.remove_vertex(graph.vertex(i))
        else:
            assert i_new >= 0
            old_to_new_nodes[i] = i_new
            i_new -= 1
    # write genome graph to file
    graph.save("genome_graph.gt")
    write_genome_gfa(graph, final_paths, old_to_new_nodes, "genome_graph.gfa")
    return


def write_genome_gfa(graph, paths, old_to_new_nodes, outfile):
    """
    Writes the final genome variation graph to a GFA file, including nodes,
    edges, and haplotype paths. Returns nothing.
    """
    with open(outfile, 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        node_lengths = {}
        for v in graph.vertices():
            f.write("S\t{}\t{}\n".format(int(v), graph.vp.seq[v]))
            node_lengths[int(v)] = len(graph.vp.seq[v])
            for w in v.out_neighbors():
                f.write("L\t{}\t+\t{}\t+\t0M\n".format(int(v), int(w)))
        for i, p in enumerate(paths):
            f.write("P\t{}\t".format(i))
            f.write("+,".join([str(old_to_new_nodes[v]) for v in p]) + "+\t")
            f.write("M,".join([str(node_lengths[old_to_new_nodes[v]]) for v in p]) + "M")
            f.write("\n")
    return

def process_output(x, diff_counter, haps, min_cov, fasta):
    """
    Analyzes Gurobi output and applies the minimal strain abundance threshold,
    thus setting any estimated abundances below this threshold to zero.
    Returns a list of final haplotype abundances.
    """
    sorted_sol = sorted(enumerate(x), key=lambda x: x[1], reverse=True)
    top_x = 20
    print("\nTop {} most abundandant strains:\n".format(top_x))
    print("\tStrain\t\tFrequency\tLength")
    count_cov = 0
    count_total = 0
    total_cov = sum([cov if cov >= min_cov else 0 for cov in x])
    final_haps_idx = []
    for idx, value in sorted_sol:
        l = len(haps[idx])
        if value >= min_cov:
            count_cov += 1
            final_haps_idx.append(idx)
            if count_total < top_x:
                print("\t{}\t\t{:.3f}\t\t{}".format(idx, value/total_cov, l))
        count_total += 1
    print("\nNumber of strains of more than {}x coverage = {}.\n".format(
            min_cov, count_cov))

    final_abundances = []
    with open(fasta, 'w') as f_out:
        l = 0
        for hap in haps:
            if l in final_haps_idx:
                cov = int(x[l])
                final_abundances.append(cov)
                f_out.write(">path{} {}x frequency={:.3f}\n".format(l, cov, cov/total_cov))
                f_out.write("{}\n".format(hap))
            else:
                final_abundances.append(0)
            l += 1

    if diff_counter != []:
        print("Node difference counts between top {} paths:\n".format(top_x))
        npaths = len(haps)
        for i in range(min(npaths, top_x)): # note: diff_counter is not sorted!
            idx_i, value_i = sorted_sol[i]
            diff_list = ["\t"]
            if x[idx_i] < min_cov:
                continue
            for j in range(min(i+1, top_x)):
                idx_j, value_j = sorted_sol[j]
                if x[idx_j] < min_cov:
                    continue
                diff_list.append("{:d}".format(int(diff_counter[i,j])))
            print("\t".join(diff_list).expandtabs(4))
        print()

    return final_abundances


def write_intermediate_files(graph, paths, adj_out):
    """Write trimmed contigs, their node paths and the overlap graph to file."""
    contig_seqs = {}
    with open("trimmed_contigs.fasta", 'w') as f:
        with open("trimmed_contigs.paths", 'w') as f2:
            for c, path in paths.items():
                seq = ""
                path_string = ""
                for v in path:
                    seq += graph.vp.seq[graph.vertex(v)]
                    path_string += "{},".format(v)
                f.write(">{}\n{}\n".format(c, seq))
                f2.write(">{}\n{}\n".format(c, path_string[:-1]))
                contig_seqs[c] = seq
    # write contig-adjacency graph to gfa
    with open("trimmed_contigs.gfa", 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        for c, adj_list in adj_out.items():
            f.write("S\t{}\t{}\n".format(c, contig_seqs[c]))
            for w in adj_list:
                f.write("L\t{}\t+\t{}\t+\t0M\n".format(c, w))


def check_containment(s1, s2):
    # check if s1 is contained in s2 while ignoring 'N' characters
    def contains(s1, s):
        idx = 0
        for seq in split1:
            old_idx = idx
            try:
                idx = s.index(seq)
                if idx < old_idx:
                    return False
                idx += len(seq)
            except ValueError:
                return False
        return True
    # build all possible s2 variations and check containment
    split1 = s1.split('N')
    split2 = s2.split('N')
    s2_variations = [""]
    for seq in split2:
        new_s2_variations = []
        for current_seq in s2_variations:
            for nuc in 'ACTG':
                extended_seq = current_seq + nuc + seq
                if contains(s1, extended_seq):
                    return True
                new_s2_variations.append(extended_seq)
        s2_variations = new_s2_variations
    return False


def build_adj_info(graph, paths, nvert, low_ab_nodes, trim_size):
    """
    Trims contig extremities and finds feasible contig concatenations.
    """
    # trim contigs
    trimmed_list = []
    delete_paths = []
    for c, path in paths.items():
        trim_left = 0
        while len(path) > 0:
            v = path[0]
            l = len(graph.vp.seq[graph.vertex(v)])
            if trim_left + l <= trim_size or v in low_ab_nodes:
                path = path[1:]
                contigs = list(graph.vp.contigs[graph.vertex(v)])
                contigs.remove(str(c))
                graph.vp.contigs[graph.vertex(v)] = contigs
                trim_left += l
            else:
                # left-trimming done
                break
        trim_right = 0
        while len(path) > 0:
            v = path[-1]
            l = len(graph.vp.seq[graph.vertex(v)])
            if trim_right + l <= trim_size or v in low_ab_nodes:
                path = path[:-1]
                contigs = list(graph.vp.contigs[graph.vertex(v)])
                contigs.remove(c)
                graph.vp.contigs[graph.vertex(v)] = contigs
                trim_right += l
            else:
                # right-trimming done
                break

        if len(path) == 0:
            delete_paths.append(c)
        elif trim_left > 0 or trim_right > 0:
            trimmed_list.append(c)
        paths[c] = path

    for c in delete_paths:
        del paths[c]

    # remove inclusions (due to trimming)
    # NOTE: very naive implementation, improve by traversing graph once
    del_set = set()
    if trim_size > 0:
        del_list = []
        print("removing inclusions due to trimming")
        path_strings = {c : ''.join([graph.vp.seq[graph.vertex(v)] for v in path])
                            for c, path in paths.items()}
        # path_strings = {c : "_".join([str(x) for x in path])
        #                     for c, path in paths.items()}
        # print(path_strings)
        # for c1 in trimmed_list:
        #     path1 = path_strings[c1]
        for c1, path1 in path_strings.items():
            for c2, path2 in path_strings.items():
                if c1 == c2:
                    continue
                elif (path1 != path2 or c1 < c2 or c2 not in trimmed_list):
                    if path1 in path2:
                        del_list.append(c1)
                    # elif check_containment(path1, path2):
                    #     del_list.append(c1)
        del_set = set(del_list)
        for c in del_set:
            del paths[c]

    starts_per_node = {v : [] for v in range(nvert)}
    for c, path in paths.items():
        v = path[0]
        starts_per_node[v].append(c)

    # build overlap list per contig (out-edges)
    dup_list = get_dup_nodes(graph)

    adj_out = {}
    start_info = {}
    for c, path in paths.items():
        adj_list, start_indexes = get_adj_out(c, paths, graph, low_ab_nodes,
                starts_per_node, dup_list, del_set)
        adj_list = [x for x in adj_list if x not in del_set]
        adj_out[c] = adj_list
        start_info[c] = {x : start_indexes[x] for x in adj_list}
    ord_adj_out = OrderedDict(sorted(adj_out.items(), key=lambda x:x[0]))
    return graph, paths, ord_adj_out, start_info


def contig_paths_to_nodes(cg_paths, g, paths, start_info, low_ab_nodes, abundances, min_ab):
    """
    Translates a list of contig paths (through the contig concatenation graph)
    to the corresponding node paths through the contig variation graph.
    Returns a list of node paths and a list of the corresponding haplotype
    sequences.
    """
    # print("start_info:")
    # print(start_info)
    print("Translating contig paths to nodes in variation graph")
    node_paths = []
    haps = []
    for contig_path in cg_paths:
        node_path = []
        idx_v_to_correct = []
        current_idx = 0
        i = 0
        c_prev = contig_path[0]
        node_path += paths[c_prev]
        for c in contig_path[1:]:
            start_idx = start_info[c_prev][c]
            current_idx += (start_idx - i)
            # traverse overlap to find low-abundance nodes
            p = paths[c]
            assert len(p) > 0
            i = 0
            while current_idx < len(node_path) and i < len(p):
                assert current_idx < len(node_path)
                if node_path[current_idx] != p[i]:
                    idx_v_to_correct.append(current_idx)
                current_idx += 1
                i += 1
            # add remaining path sequence
            node_path += p[i:]
            c_prev = c
        # reconstruct the corresponding haplotype sequence
        hap = []
        for v in node_path:
            hap.append(g.vp.seq[g.vertex(v)])
        # correct low-abundance nodes
        # corrected_path, hap = correct_path(g, idx_v_to_correct, node_path,
        #                 hap, contig_path, abundances, min_ab)
        corrected_path = node_path
        hap = ''.join(hap)
        haps.append(hap)
        node_paths.append(corrected_path)

    # nvert = len(B[0])
    # ncontigs = len(B)
    # node_paths = []
    # haps = []
    # for contig_path in cg_paths:
    #     node_path = []
    #     hap = []
    #     idx_v_to_correct = []
    #     c_prev = -1
    #     for c in contig_path:
    #         if c_prev >= 0:
    #             # find rightmost common node between paths
    #             start_node = -1
    #             for x in reversed(range(nvert)):
    #                 if B[c][x] and B[c_prev][x]:
    #                     start_node = x
    #                     break
    #             if start_node == -1:
    #                 print("ERROR: no common node found in overlap. Exiting.")
    #                 sys.exit(1)
    #         else:
    #             start_node = 0
    #
    #         # translate contigs to nodes
    #         for i, is_in_contig in enumerate(B[c]):
    #             if not is_in_contig:
    #                 continue
    #             elif i in node_path:
    #                 continue
    #             elif i < start_node:
    #                 continue
    #             elif i in low_ab_nodes:
    #                 idx_v_to_correct.append(len(node_path))
    #             node_path.append(i)
    #             hap.append(g.vp.seq[g.vertex(i)])
    #     corrected_path, hap = correct_path(g, B, idx_v_to_correct, node_path,
    #                 hap, contig_path, abundances, min_ab)
    #     node_paths.append(corrected_path)
    #     haps.append(hap)

    # filter out any duplicate paths
    print("Filtering out equivalent paths")
    remove_path_idx = []
    for h1 in range(len(haps)):
        hap1 = haps[h1]
        for h2 in range(h1+1, len(haps)):
            hap2 = haps[h2]
            if hap1 in hap2:
                remove_path_idx.append(h1)
            elif hap2 in hap1:
                remove_path_idx.append(h2)
    for idx in sorted(list(set(remove_path_idx)), reverse=True):
        node_paths.pop(idx)
        haps.pop(idx)

    print("Reduced #paths from {} to {}.".format(len(cg_paths), len(node_paths)))

    with open("haps.fasta", 'w') as f_haps:
        for i, hap in enumerate(haps):
            f_haps.write(">path{}\n{}\n".format(i, hap))
    return node_paths, haps


def enumerate_contig_paths(graph, adj_out, discarded):
    """
    Enumerates all candidate paths through the contig concatenation graph,
    allowing only feasible concatenations as specified by adj_out and avoiding
    any discarded overlaps.
    Returns a list of candidate paths.
    """
    # ncontigs = len(adj_out)
    # remove discarded overlaps from adjacency lists
    for (v1, v2) in discarded:
        adj_out[v1].remove(v2)
    # translate out-edges to in-edges
    adj_in = {c : [] for c in adj_out.keys()}
    for c, adj_list in adj_out.items():
        for i in adj_list:
            adj_in[i].append(c)
    # remove transitive edges
    remove_in = {i : [] for i in adj_out.keys()}
    for i, adj_list in adj_in.items():
        in_neighbor_count = len(adj_list)
        if in_neighbor_count < 2:
            # not enough edges for transitivity
            continue
        # check pairs of neighbors for forming a 3-clique -> transitive
        for idx1 in range(in_neighbor_count):
            v1 = adj_list[idx1]
            for idx2 in range(idx1+1, in_neighbor_count):
                v2 = adj_list[idx2]
                if v2 in adj_out[v1]:
                    remove_in[i].append(v1)
                elif v1 in adj_out[v2]:
                    remove_in[i].append(v2)
    for i, remove_list in remove_in.items():
        for v in set(remove_list):
            adj_in[i].remove(v)
            adj_out[v].remove(i)
    # find source and sink nodes
    source_nodes = []
    sink_nodes = []
    paths_to_process = []
    final_paths = []
    for v in source_nodes:
        paths_to_process.append([v])
    for c, adj_list in adj_out.items():
        if len(adj_list) > 0 and len(adj_in[c]) == 0:
            source_nodes.append(c)
            paths_to_process.append([c])
        elif len(adj_list) == 0 and len(adj_in[c]) > 0:
            sink_nodes.append(c)
        elif len(adj_list) == 0 and len(adj_in[c]) == 0:
            final_paths.append([c])
    # sink_nodes = [len(adj_out[i]) == 0 and len(adj_in[i]) > 0
    #                 for i in range(ncontigs)]
    # source_nodes = [len(adj_out[i]) > 0 and len(adj_in[i]) == 0
    #                 for i in range(ncontigs)]
    print("source contigs:", source_nodes)
    print("sink contigs:", sink_nodes)
    # traverse graph to enumerate contig paths
    while len(paths_to_process) > 0:
        path = paths_to_process.pop(0) # BFS
        v = path[-1]
        if adj_out[v] == []:
            final_paths.append(path)
        else:
            is_end = True
            for w in adj_out[v]:
                if w in path:
                    continue
                is_end = False
                paths_to_process.append(path + [w])
            if is_end:
                final_paths.append(path)
    print("{} paths found in contig graph".format(len(final_paths)))
    # print(final_paths)
    return final_paths


def check_branches2(graph, adj_out):
    """
    Check branches in overlap graph: if any branching edge is supported by a
    clique of size 3, discard all overlaps that are not supported by such a
    clique. Returns a list of contig overlaps to be discarded during path
    enumeration.
    """
    print("check_branches2")

    # translate out-edges to in-edges
    adj_in = {c : [] for c in adj_out.keys()}
    for c, adj_list in adj_out.items():
        for i in adj_list:
            adj_in[i].append(c)

    total_discard_list = []
    # for every contig with in-degree > 2, check if it is part of a 3-clique;
    for i, adj_list in adj_in.items():
        in_neighbor_count = len(adj_list)
        if in_neighbor_count <= 2:
            # not enough edges for 3-clique AND branch
            continue
        keep_list = []
        # check pairs of neighbors for forming a 3-clique
        for idx1 in range(in_neighbor_count):
            v1 = adj_list[idx1]
            for idx2 in range(idx1+1, in_neighbor_count):
                v2 = adj_list[idx2]
                if ( v2 in adj_out[v1] ) or ( v1 in adj_out[v2] ):
                    # 3-clique found
                    keep_list += [v1, v2]
            if len(adj_out[v1]) == 1:
                # don't discard out-edges for contigs with outdegree 1
                keep_list.append(v1)
        if keep_list != []:
            discarded = set(adj_list) - set(keep_list)
            for v in discarded:
                total_discard_list.append((v, i))
    return total_discard_list


def get_dup_nodes(g):
    """Returns a list of duplicate nodes (i.e. parallel nodes with identical
    sequence) in g."""
    # print("get_dup_nodes")
    dup_list = []
    for v in g.vertices():
        if v.out_degree() > 1:
            base_to_node = {}
            for w in v.out_neighbors():
                b = g.vp.seq[w]
                if b in base_to_node:
                    # dup_list.append(w)
                    base_to_node[b].append(w)
                else:
                    base_to_node[b] = [w]
            for b, node_list in base_to_node.items():
                if len(node_list) > 1:
                    dup_list += node_list
                    # print(node_list)
    dup_list = list(set(dup_list))
    print("# of duplicates: {}".format(len(dup_list)))
    return dup_list


def get_adj_out(c, paths, g, low_ab_nodes, starts_per_node, dup_nodes, del_set):
    """
    Returns an adjacency list of all feasible concatenations for a given contig.
    """
    path = paths[c]
    nvert = len(list(g.vertices()))
    start_node = path[0]
    active_contigs = set(g.vp.contigs[g.vertex(start_node)])
    rescued_contigs = set()
    start_indexes = {}
    assert len(active_contigs) > 0
    contig_seq = [g.vp.seq[g.vertex(v)] for v in path]
    # find out-neighbors
    for i, node in enumerate(path):
        old_active_contigs = active_contigs
        # update active contig set
        if node in low_ab_nodes + dup_nodes:
            assert len(path) > 1
            if i > 0:
                prev_node = path[i-1]
                alt_nodes = list(g.vertex(prev_node).out_neighbors())
            else:
                next_node = path[i+1]
                alt_nodes = list(g.vertex(next_node).in_neighbors())
            for alt in alt_nodes:
                active_contigs |= set(starts_per_node[alt])
        else:
            current_contigs = set(g.vp.contigs[g.vertex(node)])
            starts = set(starts_per_node[node])
            active_contigs = active_contigs.intersection(current_contigs)
            active_contigs |= starts
        # for any discarded active contigs, check if sequence is truly different
        for c_old in old_active_contigs - active_contigs:
            if c_old in del_set:
                continue
            idx = start_indexes[c_old]
            overlap = ''.join(contig_seq[idx:])
            assert len(overlap) > 0
            overlap_idx = 0
            keep_c = True
            overlap2 = ""
            for u in paths[c_old]:
                if not keep_c:
                    break
                u_seq = g.vp.seq[g.vertex(u)]
                overlap2 += u_seq
                for nuc in u_seq:
                    if overlap_idx >= len(overlap):
                        break
                    elif nuc == 'N' or overlap[overlap_idx] == 'N':
                        overlap_idx += 1
                        continue
                    elif nuc != overlap[overlap_idx]:
                        keep_c = False
                        break
                    overlap_idx += 1
            if keep_c:
                rescued_contigs.add(c_old)
                # print("rescued contig overlap {} -> {}".format(c, c_old))
        # keep track of start points
        for x in active_contigs:
            if x not in start_indexes:
                start_indexes[x] = i
    assert c in active_contigs
    active_contigs |= rescued_contigs
    active_contigs.remove(c)
    return list(active_contigs), start_indexes

    # start_node = path.index(1)
    # active_contigs = set(int(x) for x in g.vp.contigs[g.vertex(start_node)])
    # assert len(active_contigs) > 0
    # prev_node = -1
    # for node in range(nvert):
    #     if path[node] == 0:
    #         # node not in path
    #         continue
    #     elif node in low_ab_nodes + dup_nodes:
    #         if prev_node >= 0:
    #             alt_nodes = list(g.vertex(prev_node).out_neighbors())
    #         else:
    #             next_node = node + 1 + path[node+1:].index(1)
    #             alt_nodes = list(g.vertex(next_node).in_neighbors())
    #         for alt in alt_nodes:
    #             active_contigs |= set(starts_per_node[alt])
    #     else:
    #         current_contigs = set(int(x) for x in g.vp.contigs[g.vertex(node)])
    #         starts = set(starts_per_node[node])
    #         active_contigs = active_contigs.intersection(current_contigs)
    #         active_contigs |= starts
    #     prev_node = node
    # assert c in active_contigs
    # active_contigs.remove(c)
    # return list(active_contigs)


def get_extremity_nodes(graph, min_depth):
    """
    Returns a list of node IDs corresponding to extremity nodes in the input
    graph.
    """
    # list all nodes less than <min_depth> from the source/sink nodes
    #graph = load_graph(graphfile)
    path_nodes = []
    if min_depth == 0:
        return path_nodes
    # traverse graph starting from extremities
    for v in graph.vertices():
        if v.in_degree() == 0:
            # forward DFS until depth > min_depth
            forward = True
            path = []
            depth = 0
            [path, depth] = dfs(graph, v, min_depth, depth, path, forward)
            path_nodes += path
        if v.out_degree() == 0:
            # backward DFS until depth > min_depth
            forward = False
            path = []
            depth = 0
            [path, depth]= dfs(graph, v, min_depth, depth, path, forward)
            path_nodes += path
    return list(set(path_nodes))


def dfs(graph, v, min_depth, depth, path, forward):
    """
    DFS for extremity nodes. Returns the list of nodes traversed until
    min_depth was reached and the depth at which the search was terminated.
    """
    depth += len(graph.vp.seq[v])
    if depth >= min_depth:
        return [path, depth]
    path.append(int(v))
    # continue DFS
    if forward:
        for w in v.out_neighbors():
            [path, depth] = dfs(graph, w, min_depth, depth, path, forward)
    else:
        for w in v.in_neighbors():
            [path, depth] = dfs(graph, w, min_depth, depth, path, forward)
    return [path, depth]


def correct_path(g, idx_v_to_correct, path, hap, contigs, abundances, min_ab):
    """
    Replaces all low-abundance nodes in the input path by their best
    alternatives.
    Returns the corrected path and the corresponding haplotype sequence.
    """
    # nvert = len(B[0])
    # ncontigs = len(B)
    for i in range(len(idx_v_to_correct)):
        idx = idx_v_to_correct[i]
        v = path[idx]
        if idx == 0:
            v_next = path[idx+1]
            candidate_nodes = list(g.vertex(v_next).in_neighbors())
        elif idx == len(path)-1 or ( i+1<len(idx_v_to_correct) and
                                     idx_v_to_correct[i+1] == idx+1 ):
            v_prev = path[idx-1]
            candidate_nodes = list(g.vertex(v_prev).out_neighbors())
        else:
            v_next = path[idx+1]
            v_prev = path[idx-1]
            candidate_nodes = list(
                                set(g.vertex(v_prev).out_neighbors()) &
                                set(g.vertex(v_next).in_neighbors())
            )
        max_ab = abundances[int(v)]
        best_node = v
        best_seq = g.vp.seq[v]
        for node in candidate_nodes:
            if abundances[int(node)] <= max_ab:
                continue
            node_seq = g.vp.seq[node]
            node_ab = abundances[int(node)]
            if (node_seq == best_seq
                    and node_ab + abundances[int(v)] >= min_ab):
                # equivalent node, allow replacement independent of contigs
                max_ab = node_ab
                best_node = node
                best_seq = node_seq
                break
            else:
                for k in g.vp.contigs[node]:
                    if k in contigs:
                        max_ab = abundances[int(node)]
                        best_node = node
                        break
        # replace the node with its best alternative
        path[idx] = int(best_node)
        hap[idx] = best_seq
    return path, ''.join(hap)


def optimize2(a, nvert, paths, skip_nodes, max_strains, node_lengths,
                        min_cov, min_cov_final, threads):
    m = Model('qp')
    obj = QuadExpr()
    npaths = len(paths)
    P = zeros((nvert,npaths))

    x = m.addVars(list(range(npaths)), lb=0, ub=1.05*max(a), vtype=GRB.CONTINUOUS, name='x')
    X = [x[i] for i in range(npaths)]
    X = array([X]).reshape(npaths,1)

    # add indicator variables to count strains
    if max_strains > 0:
        # add indicator variables for counting strains
        x_binary = m.addVars(list(range(npaths)), vtype=GRB.BINARY, name='strain_indicator')
        for i in range(npaths):
            max_x = 2*max(a)
            m.addConstr(x_binary[i] >= (x[i]-min_cov)/max_x)
        # define total strain count
        sum_x_binary = LinExpr()
        for i in range(npaths):
            sum_x_binary += x_binary[i]
        # bound the number of strains
        m.addConstr(sum_x_binary <= max_strains)

    # Store paths in P: p_ij = 1 if node i contains path j
    print('\nSave for every node which paths are passing through:')
    for i in tqdm(range(npaths)):
        for v in paths[i]:
            P[int(v),i] = 1
    npaths = len(paths)
    del paths
    m.update()

    if max_strains > 0:
        sum_x_binary = LinExpr()
        for i in range(npaths):
            sum_x_binary += x_binary[i]

    # Store node abundance values per path combination
    path_combinations_seen = {}
    for v in tqdm(range(nvert)):
        if v in skip_nodes:
            continue
        # check if v adds a new combination of paths and store abundance
        path_combi = tuple(P[v])
        try:
            path_combinations_seen[path_combi].append([a[v], node_lengths[v]])
        except KeyError:
            path_combinations_seen[path_combi] = [[a[v], node_lengths[v]]]

    n_combi = len(path_combinations_seen)
    y = m.addVars(list(range(n_combi)), lb=0, vtype=GRB.CONTINUOUS, name='y')

    # Define the objective
    print('\nDefine the objective function:')
    n_eval = 0
    combi_idx = -1
    for path_combi, ab_list in path_combinations_seen.items():
        combi_idx += 1
        # calculate the average (weighted) abundance value for this combi
        abundance = sum([x*y for [x,y] in ab_list])/sum([y for [x,y] in ab_list])
        # sum the calculated abundances of strains
        sum_xv = 0
        for idx, i in enumerate(path_combi):
            if i == 1:
                sum_xv += x[idx]
        # compute squared error
        # obj += (abundance - sum_xv) * (abundance - sum_xv)
        # obj += (abundance - sum_xv) * (abundance - sum_xv) / max(0.00001, abundance)

        # # absolute difference
        obj += y[combi_idx] # abs(abundance - sum_xv)
        # obj += y[combi_idx] / max(0.00001, abundance)
        # set constraints on y[v] to obtain absolute value
        m.addConstr(y[combi_idx] >= sum_xv - abundance, "y_{}_-".format(combi_idx))
        m.addConstr(y[combi_idx] >= abundance - sum_xv, "y_{}_+".format(combi_idx))
        n_eval += 1
    assert n_eval > 0
    obj *= (1/n_eval)

    # set objective and minimize
    m.setObjective(obj, GRB.MINIMIZE)
    print('\nObjective function ready, starting Gurobi optimization:\n')
    m.update()

    m.Params.LogToConsole = 1
    m.Params.Threads = threads
    m.Params.NumericFocus = 0
    m.Params.PoolSearchMode = 0
    m.Params.PoolSolutions = 10
    m.Params.Method = 4

    #Minimize the model for the given objective function and constraints
    print("\n*** Phase 1 optimization***\n")
    m.optimize()

    print("Number of solutions = {}".format(m.solcount))
    if m.status == GRB.Status.OPTIMAL:
        x_sol = []
        nstrains_test = 0
        for v in m.getVars():
            if 'x' in v.varName:
                x_sol.append(v.x)
        print('\nObjective value: %g' % m.objVal)
        objVal = m.objVal

        selected_strains = [1 if cov > min_cov_final else 0 for cov in x_sol]
        nstrains = sum(selected_strains)
        print("#strains / #paths = {} / {}".format(nstrains, npaths))

        # run phase 2 optimization:
        print("\n*** Phase 2 optimization***\n")
        m.reset()
        for i in range(npaths):
            if selected_strains[i] == 0:
                m.addConstr(x[i] == 0)
        m.optimize()

        print("\n***Final optimization results:***")
        x_final = []
        nstrains_test = 0
        for v in m.getVars():
            if 'x' in v.varName:
                x_final.append(v.x)
        print('\nObjective value: %g' % m.objVal)
        objVal = m.objVal

        selected_strains = [1 if cov > min_cov_final else 0 for cov in x_final]
        nstrains = sum(selected_strains)
        print("#strains / #paths = {} / {}".format(nstrains, npaths))

        #Comparison of the paths, see how many nodes are different between paths
        if npaths > 5:
            diff_counter = []
        else:
            diff_counter = zeros((npaths,npaths))
            for i in tqdm(range(npaths)):
                for j in range(npaths):
                    for k in range(nvert):
                        if P[k,i] != P[k,j]:
                            diff_counter[i,j] += 1

        return(x_final, objVal, diff_counter)

    else:
        try:
            m.computeIIS()
            # Print the names of all of the constraints in the IIS set.
            print("IIS constraints:")
            for c in m.GetConstrs():
                if c.Get(GRB.IntAttr.IISConstr) > 0:
                    print(c.Get(GRB.StringAttr.ConstrName))
            # Print the names of all of the variables in the IIS set.
            print("IIS variables:")
            for v in m.GetVars():
                if v.Get(GRB.IntAttr.IISLB) > 0 or v.Get(GRB.IntAttr.IISUB) > 0:
                    print(v.Get(GRB.StringAttr.VarName))
            print("ERROR: Infeasible model.")
        except:
            #print(m.getAttr(GRB.Attr.UnbdRay, m.getVars()))
            print("ERROR: Unbounded model.")
        print('\nNo optimal solution found, exiting.')
        sys.exit(1)


def optimize(a, nvert, paths, skip_nodes, reduce_obj, max_strains,
                        min_cov, min_cov_final, threads):
    """
    Defines Gurobi minimization problem and then applies the LP solver.
    Returns the solution values, the corresponding objective value, and a matrix
    counting pairwise differences between haplotypes.
    """
    t1_start = time.perf_counter()
    t2_start = time.process_time()
    obj_func = 1
    #Define a new model
    if obj_func < 4:
        m = Model('lp')
        obj = LinExpr()
    else:
        m = Model('qp')
        obj = QuadExpr()

    npaths = len(paths)          #The number of feasible paths
    P = zeros((nvert,npaths))   #Matrix determining which node is in which path

    x = m.addVars(list(range(npaths)), lb=0, ub=1.05*max(a), vtype=GRB.CONTINUOUS, name='x')
    X = [x[i] for i in range(npaths)]
    X = array([X]).reshape(npaths,1)    #Set x in an array for multiplication

    # If objective involves absolute values, add extra variables
    if obj_func < 4 or obj_func == 6:
        y = m.addVars(list(range(nvert)), lb=0, vtype=GRB.CONTINUOUS, name='y')
    # add indicator variables to count strains
    if max_strains > 0:
        # add indicator variables for counting strains
        x_binary = m.addVars(list(range(npaths)), vtype=GRB.BINARY, name='strain_indicator')
        for i in range(npaths):
            max_x = 2*max(a)
            m.addConstr(x_binary[i] >= (x[i]-min_cov)/max_x)
        # define total strain count
        sum_x_binary = LinExpr()
        for i in range(npaths):
            sum_x_binary += x_binary[i]
        # bound the number of strains
        m.addConstr(sum_x_binary <= max_strains)

    # Store paths in P: p_ij = 1 if node i contains path j
    print('\nSave for every node which paths are passing through:')
    for i in tqdm(range(npaths)):
        for v in paths[i]:
            P[int(v),i] = 1
    npaths = len(paths)
    del paths
    m.update()

    if max_strains > 0:
        sum_x_binary = LinExpr()
        for i in range(npaths):
            sum_x_binary += x_binary[i]

    # Define the objective function
    path_combinations_seen = {}
    count_reduce_effect = 0
    print('\nDefine the objective function:')
    n_eval = 0
    for v in tqdm(range(nvert)):
        if v in skip_nodes:
            if obj_func < 4 or obj_func == 6:
                m.addConstr(y[v] == 0, "y_{}_0".format(v))
            continue
        if reduce_obj > 0:
            # check if v adds a new combination of paths
            path_combi = tuple(P[v])
            try:
                combi_count = path_combinations_seen[path_combi]
                if combi_count >= reduce_obj:
                    count_reduce_effect += 1
                    continue
                else:
                    path_combinations_seen[path_combi] += 1
            except KeyError:
                path_combinations_seen[path_combi] = 1
        # sum the calculated abundances of strains through v
        #sum_xv = dot(P[v,:],X)[0] # memory expensive!!!
        sum_xv = 0
        for idx, i in enumerate(P[v]):
            if i == 1:
                sum_xv += x[idx]
        #abundance = int(round(a[v]))
        abundance = a[v]
        if abundance == 0:
            continue
        if obj_func == 1 or obj_func == 6:
            # absolute difference
            obj += y[v] #abs(abundance - sum_xv)
            # set constraints on y[v] to obtain absolute value
            m.addConstr(y[v] >= sum_xv - abundance, "y_{}_-".format(v))
            m.addConstr(y[v] >= abundance - sum_xv, "y_{}_+".format(v))
        elif obj_func == 2:
            # relative absolute difference (linear)
            obj += y[v] / abundance # abs(abundance - sum_xv) / abundance
            # set constraints on y[v] to obtain absolute value
            m.addConstr(y[v] >= sum_xv - abundance, "y_{}_-".format(v))
            m.addConstr(y[v] >= abundance - sum_xv, "y_{}_+".format(v))
        elif obj_func == 3:
            # relative absolute difference (sqrt)
            obj += y[v] / math.sqrt(abundance) # abs(abundance - sum_xv) / math.sqrt(abundance)
            # set constraints on y[v] to obtain absolute value
            m.addConstr(y[v] >= sum_xv - abundance, "y_{}_-".format(v))
            m.addConstr(y[v] >= abundance - sum_xv, "y_{}_+".format(v))
        elif obj_func == 4 or obj_func == 7:
            # squared error
            obj += (abundance - sum_xv) * (abundance - sum_xv)
        elif obj_func == 5:
            # relative squared error
            obj += (abundance - sum_xv) * (abundance - sum_xv) / max(0.00001, abundance)
        else:
            print("Objective function not recognized; exiting.")
            sys.exit(1)
        n_eval += 1
    assert n_eval > 0

    if obj_func == 6 or obj_func == 7:
        # add strain penalty
        sum_x_binary = LinExpr()
        for i in range(npaths):
            sum_x_binary += x_binary[i]
        if obj_func == 6:
            obj *= (1/n_eval)*(npaths+sum_x_binary)
        else:
            obj += sum_x_binary*math.log(n_eval)
    else:
        obj *= (1/n_eval)

    # set objective and minimize
    m.setObjective(obj, GRB.MINIMIZE)
    print('\nObjective function ready, starting Gurobi optimization:\n')
    if reduce_obj > 0:
        print("Reduced objective function by {} nodes (nvert = {}).\n".format(count_reduce_effect, nvert))
    m.update()

    m.Params.LogToConsole = 1
    m.Params.Threads = threads
    m.Params.NumericFocus = 0
    m.Params.PoolSearchMode = 0
    m.Params.PoolSolutions = 10
    m.Params.Method = 4
    #m.Params.InfUnbdInfo = 1
    #m.Params.PreQLinearize = 1 # default = automatic

    #Minimize the model for the given objective function and constraints
    print("\n*** Phase 1 optimization***\n")
    m.optimize()

    print("Number of solutions = {}".format(m.solcount))
    if m.status == GRB.Status.OPTIMAL:
        x_sol = []
        nstrains_test = 0
        for v in m.getVars():
            if 'x' in v.varName:
                x_sol.append(v.x)
                # print("strain abundance: ", v.x)
            elif obj_func >= 6 and 'strain_indicator' in v.varName:
                nstrains_test += int(v.x)
                # print("strain indicator: ", v.x)
        print('\nObjective value: %g' % m.objVal)
        objVal = m.objVal

        selected_strains = [1 if cov > min_cov_final else 0 for cov in x_sol]
        nstrains = sum(selected_strains)
        print("#strains / #paths = {} / {}".format(nstrains, npaths))
        #print("nstrains_test = {}".format(nstrains_test))
        #assert nstrains == nstrains_test

        # run phase 2 optimization:
        print("\n*** Phase 2 optimization***\n")
        m.reset()
        for i in range(npaths):
            if selected_strains[i] == 0:
                m.addConstr(x[i] == 0)
        m.optimize()

        print("\n***Final optimization results:***")
        x_final = []
        nstrains_test = 0
        for v in m.getVars():
            if 'x' in v.varName:
                x_final.append(v.x)
                # print("strain abundance: ", v.x)
            elif obj_func >= 6 and 'strain_indicator' in v.varName:
                nstrains_test += int(v.x)
                # print("strain indicator: ", v.x)
        print('\nObjective value: %g' % m.objVal)
        objVal = m.objVal

        selected_strains = [1 if cov > min_cov_final else 0 for cov in x_final]
        nstrains = sum(selected_strains)
        print("#strains / #paths = {} / {}".format(nstrains, npaths))

        # with open('strain_abundance.txt', 'w') as f:
        #     for i in range(len(x)):
        #         f.write('strain_%s:%g\n'%(i,x[i]))

        #Comparison of the paths, see how many nodes are different between paths
        if npaths > 5:
            diff_counter = []
        else:
            diff_counter = zeros((npaths,npaths))
            for i in tqdm(range(npaths)):
                for j in range(npaths):
                    for k in range(nvert):
                        if P[k,i] != P[k,j]:
                            diff_counter[i,j] += 1

        t1_stop = time.perf_counter()
        t2_stop = time.process_time()
        print("\nPath selection completed")
        print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
        print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
        return(x_final, objVal, diff_counter)

    else:
        try:
            m.computeIIS()
            # Print the names of all of the constraints in the IIS set.
            print("IIS constraints:")
            for c in m.GetConstrs():
                if c.Get(GRB.IntAttr.IISConstr) > 0:
                    print(c.Get(GRB.StringAttr.ConstrName))
            # Print the names of all of the variables in the IIS set.
            print("IIS variables:")
            for v in m.GetVars():
                if v.Get(GRB.IntAttr.IISLB) > 0 or v.Get(GRB.IntAttr.IISUB) > 0:
                    print(v.Get(GRB.StringAttr.VarName))
            print("ERROR: Infeasible model.")
        except:
            #print(m.getAttr(GRB.Attr.UnbdRay, m.getVars()))
            print("ERROR: Unbounded model.")

        print('\nNo optimal solution found, exiting.')
        sys.exit(1)


def cyclic(adj_out):
    """Return True if the directed input graph has a cycle."""
    graph = Graph(directed=True)
    vprop = graph.new_vertex_property('string')
    graph.vp.contig = vprop
    contig_to_node = {}
    node_to_contig = {}
    for c, adj_list in adj_out.items():
        v = graph.add_vertex()
        # print(int(v), c)
        contig_to_node[c] = v
        node_to_contig[v] = c
        graph.vp.contig[v] = c
    for c, adj_list in adj_out.items():
        for c2 in adj_list:
            graph.add_edge(contig_to_node[c], contig_to_node[c2])
    if is_DAG(graph):
        return []
    else:
        # print("CYCLES:")
        cycles = []
        for cycle in all_circuits(graph):
            cycles.append([node_to_contig[v] for v in cycle])
        return cycles


if __name__ == '__main__':
    sys.exit(main())

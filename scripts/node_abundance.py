#!/usr/bin/env python3

import sys
import json
from graph_tool.all import *
import subprocess
import argparse
from numpy import array, zeros
from tqdm import tqdm # progress tracker


__author__ = "Jasmijn Baaijens"
__license__ = "MIT"

usage = "Compute node abundances by aligning original reads to variation graph."


def main():
    parser = argparse.ArgumentParser(prog='node_abundance.py', description=usage)
    parser.add_argument('-g', '--graph', dest='graph', type=str, required=True, help='Input variation graph in graph-tool format')
    parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads in fastq format')
    parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads in fastq format')
    parser.add_argument('-c', '--contigs', dest='contigs', type=str, required=True, help='Input contigs in fastq format')
    parser.add_argument('-vg', '--vg_path', dest='vg', type=str, required=True, help="Path to vg executable")
    parser.add_argument('-t', '--threads', default=1, help="Set number of threads used for vg.")
    parser.add_argument('--skip_alignments', action='store_true', help="Read vg alignments from existing file.")
    parser.add_argument('--abundance_file', type=str, help="Read pre-computed node abundances from this file.")
    parser.add_argument('--ref', default='', help="Optional: reference-guided contig graph construction. Provide a reference genome for vg construct")
    parser.add_argument('--vcf', default='', help="Optional: reference-guided contig graph construction. Provide a VCF file for vg construct")
    args = parser.parse_args()

    graph_name = args.graph.rstrip('.gt')
    reads_forward = args.forward
    reads_reverse = args.reverse
    aln_name = "aln"
    vg = args.vg

    try:
        subprocess.check_call(vg + " version > /dev/null 2>&1", shell=True)
    except:
        print("""
ERROR: vg executable not found. Please check if the correct path was
specified in the -vg argument, and check if your vg installation works.
Exiting.
        """)
        sys.exit(1)

    if not args.skip_alignments:
        # run vg to get read alignments
        vg_map_reads(graph_name, reads_forward, reads_reverse, aln_name, vg,
            args.threads, args.ref, args.vcf, args.contigs)

    aln_file = "{}.filtered.json".format(aln_name)

    if args.ref != '' and args.vcf != '':
        graph = gfa_to_gt(graph_name, args.contigs, vg, args.threads)
    else:
        graph = load_graph(args.graph)

    if not args.abundance_file:
        # process alignments and compute abundance rates
        node_abundances, connection_counts = get_node_abundances(graph, aln_file)
        #node_abundances, edge_abundances = get_abundances(graph, aln_file)
    else:
        node_abundances = {}
        with open(args.abundance_file, 'r') as f:
            for line in f:
                [ID, ab] = line.rstrip().split(':')
                node_abundances[int(ID)] = float(ab)
        connection_counts = {}
        print("Using node abundance from {}".format(args.abundance_file))

    # store nodes on graph extremities such that these can be ignored during
    # objective evaluation; also store nodes to be removed entirely.
    low_ab_nodes = []
    abundance_file = open('node_abundance.txt', 'w')
    node_id = 0
    for node, ab in node_abundances.items():
        abundance_file.write(str(node_id) + ':' + str(ab) + '\n')
        node_id += 1
    abundance_file.close()
    return


def gfa_to_gt(graph_name, contigs, vg, threads):
    """
    Reads a graph from a GFA-file and converts it to graph-tool format.
    Returns graph in gt-format.
    """
    # Define a graph with its vertex properties
    g = Graph(directed=True)
    vprop = g.new_vertex_property('string')
    g.vp.seq = vprop
    vprop = g.new_vertex_property('vector<string>')
    g.vp.contigs = vprop

    # read gfa and add vertices and edges to graph
    with open(graph_name + ".gfa", 'r') as f:
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
            elif line[0] == 'L':
                # edge
                assert line[2] == line[4] == "+"
                v1 = int(line[1]) - 1
                v2 = int(line[3]) - 1
                g.add_edge(v1, v2)

    # map contigs back to graph to find their paths
    contig_aln = "contigs_aln"
    kmer_size = 11
    subprocess.check_call( vg +
        " map -f {0} -f {0} -x {1}.xg -g {1}.gcsa -k {2} -t {3} > {4}.gam".format(
        contigs, graph_name, kmer_size, threads, contig_aln),
        shell=True)
    # # filter alignments
    # subprocess.check_call( vg +
    #     " filter -r 0.90 -s 2 -fu -t {1} {0}.gam > {0}.filtered.gam".format(
    #     contig_aln, threads),
    #     shell=True)
    # output json
    print("convert gam to json format")
    subprocess.check_call( vg +
        " view -a {0}.gam > {0}.json".format(contig_aln),
        shell=True)

    # subprocess.check_call("rm -f {}.json".format(contig_aln), shell=True)
    # contig_dict = {}
    # with open(contigs, 'r') as f:
    #     c = 0
    #     for line in f:
    #         if line[0] == '>':
    #             continue
    #         seq = line.rstrip()
    #         aln_line = subprocess.check_call( vg +
    #                 " align -Q {} -s {} {}.vg >> {}.gam".format(
    #                 c, seq, graph_name, contig_aln), shell=True)
    #         contig_dict[c] = seq
    #         c += 1

    # read contigs and trace path through graph
    with open(contig_aln + '.json', 'r') as f:
        c = 0
        seq_id_prev = ""
        for line in f:
            aln = json.loads(line)
            seq_id, path = get_contig_aln_path(aln)
            if seq_id == seq_id_prev:
                continue
            seq_id_prev = seq_id
            print(c)
            print(path)
            v_old = -1
            for v in path:
                g.vp.contigs[g.vertex(v)].append(str(c))
                if v_old >= 0 and v not in g.vertex(v_old).out_neighbors():
                    print("Adding extra edge!")
                    g.add_edge(v_old, v)
                v_old = v
            c += 1

    # save graph in gt format
    g.save(graph_name + '.gt')
    return g


def get_contig_aln_path(aln):
    """Returns the node-path specified by a given GAM (json) alignment."""
    seq_id = aln["name"]
    seq = aln["sequence"]
    try:
        path = aln["path"]
        mapping = path["mapping"]
    except KeyError:
        # unmapped
        print("WARNING: contig unmapped!!!")
        return seq_id, []
    node_path = []
    for node_info in mapping:
        position = node_info["position"]
        node_id = position["node_id"]
        try:
            is_reverse = position["is_reverse"]
        except KeyError:
            is_reverse = False
        node = int(node_id)-1
        node_path.append(node)
    if is_reverse:
        return seq_id, list(reversed(node_path))
    else:
        return seq_id, node_path


def vg_map_reads(graph_name, reads_forward, reads_reverse, aln_name, vg,
        threads, ref, vcf, contigs):
    """
    Executes all vg steps to map the given paired-end reads to the contig
    variation graph and stores the alignments in a json file. Returns nothing.
    """
    kmer_size = 11
    subprocess.check_call("mkdir -p tmp", shell=True)
    # build graph in vg format
    print("building graph in vg format...")
    if ref == '' or vcf == '':
        subprocess.check_call( vg +
            " view -F -v {0}.gfa > {0}.vg".format(graph_name), shell=True)
    else:
        subprocess.check_call( vg +
            " construct -r {} -v {} > {}.tmp.vg".format(ref, vcf, graph_name), shell=True)
        # # map contigs
        # contig_aln = "contigs_aln"
        # subprocess.check_call("rm -f {}.json".format(contig_aln), shell=True)
        # with open(contigs, 'r') as f:
        #     c = 0
        #     for line in f:
        #         if line[0] == '>':
        #             continue
        #         seq = line.rstrip()
        #         aln_line = subprocess.check_call( vg +
        #                 " align -Q {} -s {} {}.tmp.vg >> {}.gam".format(
        #                 c, seq, graph_name, contig_aln), shell=True)
        #         c += 1
        # # integrate contig paths
        # subprocess.check_call( vg +
        #     " mod -i {0}.gam {1}.tmp.vg > {1}.vg".format(
        #     contig_aln, graph_name), shell=True)
        # convert graph to gfa
        subprocess.check_call( vg +
            " view {0}.vg > {0}.gfa".format(graph_name), shell=True)

    # map reads (global alignment)
    print("indexing graph...")
    #input_check = input("Is your current $TMPDIR large enough to store intermediate gcsa files?")
    subprocess.check_call( "export TMPDIR=tmp && " + vg +
        " index -x {0}.xg -g {0}.gcsa -k {1} -X 2 -Z 500 -t {2}".format(
            graph_name, kmer_size, threads) +
        " --path-only" # only index paths defined by contigs
        " {0}.vg".format(graph_name), shell=True)
    print("mapping reads to graph...")
    subprocess.check_call( vg +
        " map -f {0} -f {1} -x {2}.xg -g {2}.gcsa -k {4} -t {5} > {3}.gam".format(
        reads_forward, reads_reverse, graph_name, aln_name, kmer_size, threads),
        shell=True)

    # filter for primary alignments
    print("filter alignments...")
    subprocess.check_call( vg +
        " filter -r 0.90 -s 2 -fu -t {1} {0}.gam > {0}.filtered.gam".format(
        aln_name, threads),
        shell=True)

    # output json
    print("convert gam to json format")
    subprocess.check_call( vg +
        " view -a {0}.filtered.gam > {0}.filtered.json".format(aln_name),
        shell=True)
    return


def get_node_abundances(graph, aln_file):
    """
    Reads vg alignments and computes node abundance values and edge abundances
    (connection counts per edge).
    Returns node abundance list and connection count dict.
    """
    nodes = {}
    for node in graph.vertices():
        # update node ID, because alignments are to GFA graph which has 1-based
        # node IDs
        node_id = int(node)+1
        seq = graph.vp.seq[node]
        nodes[node_id] = seq

    bases_per_node = {} # map node IDs to read alignments
    for node in nodes:
        bases_per_node[node] = 0

    # count node connections used by read alignments
    connection_counts = {}

    #node_to_seq = {node : [] for node in nodes}

    print("Processing alignments...")
    with open(aln_file, 'r') as aln_json:
        for line in tqdm(aln_json):
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
                node2 = int(node_id)-1
                node_seq = nodes[node_id]
                node_len = len(node_seq)
                aln_len = 0
                try:
                    offset = int(position["offset"])
                except KeyError:
                    offset = 0
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
                if node0 != "":
                    # store 2-edge connection
                    if (node0, node2) in connection_counts:
                        try:
                            connection_counts[(node0, node2)][node1] += 1
                        except KeyError:
                            connection_counts[(node0, node2)][node1] = 1
                    else:
                        connection_counts[(node0, node2)] = {node1 : 1}
                if node1 != "":
                    # store direct connection
                    if (node1, node2) in connection_counts:
                        try:
                            connection_counts[(node1, node2)][-1] += 1
                        except KeyError:
                            connection_counts[(node1, node2)][-1] = 1
                    else:
                        connection_counts[(node1, node2)] = {-1 : 1}
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
        node_abundance_list[node-1] = node_abundance
        # print(node, ":", node_abundance)
        # f.write("{0}:{1}\n".format(node-1, node_abundance))

    return node_abundance_list, connection_counts


def count_contigs(contig_file):
    """Returns number of contigs in contig file."""
    count = 0
    with open(contig_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                count += 1
    return count


def trim_graph(graph, low_ab_nodes, connection_counts, contig_file):
    """
    Trim low abundance nodes from input variation graph.
    WARNING: may break contig paths. CURRENTLY NOT IN USE.
    Returns the updated graph and a list of affected contig information.
    """
    # build matrix B (nodes x contigs) indicating which node is on which contigs
    ncontigs = count_contigs(contig_file)
    nvert = len(list(graph.vertices()))
    B = zeros((nvert,ncontigs))
    for v in graph.vertices():
        for k in range(ncontigs):
            if str(k) in list(graph.vp.contigs[v]):
                B[int(v),k]=1
    # clear contigs info such that we can store it anew
    for v in graph.vertices():
        graph.vp.contigs[v] = []
    # trim graph and update contigs
    trimmed_contigs = []
    new_id = 0
    for k in range(ncontigs):
        old_id = k
        nodes = [B[i,k] for i in range(nvert)]
        seq = ""
        path = []
        identical_contigs = set()
        node0 = ""
        node1 = ""
        for i in range(nvert):
            if nodes[i] == 0:
                continue
            elif i not in low_ab_nodes: # node_abundance_list[i] >= min_abundance:
                node0 = node1
                node1 = i
                seq += graph.vp.seq[graph.vertex(i)]
                path.append(i)
                contigs = [x for x, e in enumerate(B[i]) if e != 0]
                if len(identical_contigs) > 0:
                    identical_contigs.intersection_update(contigs)
                else:
                    identical_contigs = set(contigs)
                continue
            # find the neighboring node pair (node0, next_node)
            node0 = node1
            node1 = ""
            if node0 == "":
                # erroneous node is first node of path; skip node
                continue
            elif sum(nodes[i+1:]) == 0:
                # erroneous node is last node of path; skip node
                break
            next_node = i + 1 + nodes[i+1:].index(1)
            # select best alternative node (majority vote)
            try:
                counts_dict = connection_counts[(node0, next_node)]
            except KeyError:
                counts_dict = {}
            sorted_keys = [(x, counts_dict[x]) for x in sorted(
                    counts_dict, key=counts_dict.get, reverse=True)]
            for key_value_pair in sorted_keys:
                key = key_value_pair[0]
                if key != i and key not in low_ab_nodes: # node_abundance_list[key] >= min_abundance:
                    node1 = key
                    break
            if node1 == "" or next_node in low_ab_nodes: # node_abundance_list[next_node] < min_abundance:
                # no alternative node, break path here
                print("!!! breaking contig {} at node {} !!!".format(k, i))
                if len(identical_contigs) == 1:
                    # not an inclusion so keep sequence
                    trimmed_contigs.append([new_id, seq, path])
                    # add contigs info to all nodes in path
                    for node in path:
                        v = graph.vertex(node)
                        graph.vp.contigs[v].append(str(new_id))
                    new_id += 1
                else:
                    print("inclusion contig...")
                seq = ""
                path = []
                identical_contigs = set()
                continue
            # update path and contig seq
            # print("corrected node {} by node {} in contig {}".format(i, node1, k))
            seq += graph.vp.seq[graph.vertex(node1)]
            path.append(node1)
            # update identical contig set
            contigs = [x for x, e in enumerate(B[node1]) if e != 0]
            if len(identical_contigs) > 0:
                identical_contigs.intersection_update(contigs)
            else:
                identical_contigs = set(contigs)
            continue
        # store final contig sequence
        if len(identical_contigs) == 1:
            trimmed_contigs.append([new_id, seq, path])
            # add contigs info to all nodes in path
            for node in path:
                v = graph.vertex(node)
                graph.vp.contigs[v].append(str(new_id))
            new_id += 1
    # now remove nodes with too little abundance
    for v in sorted([int(u) for u in graph.vertices()], reverse=True):
        if int(v) in low_ab_nodes: # node_abundance_list[int(v)] < min_abundance:
            graph.remove_vertex(v)
    return graph, trimmed_contigs


def get_extremity_nodes(graph, min_depth):
    """
    Returns a list of node IDs corresponding to extremity nodes in the input
    graph.
    """
    # list all nodes less than <min_depth> from the source/sink nodes
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


if __name__ == '__main__':
    sys.exit(main())

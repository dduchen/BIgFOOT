#!/usr/bin/env python3

import sys, os
import json
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort
from collections import OrderedDict
from Bio import pairwise2


def get_seq(g, v, ori):
    """Extract nucleotide sequence from graph for given node and orientation"""
    seq = g.vp.seq[g.vertex(int(v))]
    if ori == "-":
        seq = rev_comp(seq)
    return seq


def rev_comp(seq):
    """Build reverse complementary DNA sequence"""
    c_dict = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N'}
    try:
        c_seq = ''.join(c_dict[x] for x in reversed(seq.upper()))
    except KeyError as e:
        print("ERROR: segment has non-ACTGN characters")
        sys.exit(1)
    return c_seq


def has_simple_paths(g):
    """
    Check if input graph has any simple paths.
    """
    bad_edges = []
    for e in g.edges():
        if e.source().out_degree() == 1 and e.target().in_degree() == 1:
            bad_edges.append([e.source(), e.target()])
    if bad_edges != []:
        print(bad_edges)
        return True
    return False


def reorder_graph(g, paths, vertex_ordering):
    """
    Rename vertices to integers matching a given vertex ordering and
    update edges and paths accordingly.
    """
    new_g = Graph(directed=True)
    new_g.vp.seq = new_g.new_vertex_property('string')
    new_g.vp.contigs = new_g.new_vertex_property('vector<string>')
    new_g.ep.ori = new_g.new_edge_property('string')
    new_g.vp.ab = new_g.new_vertex_property('float')
    old_to_new = {}
    for i, v in enumerate(vertex_ordering):
        old_to_new[v] = i
        node = new_g.add_vertex()
        new_g.vp.seq[node] = g.vp.seq[g.vertex(v)]
        new_g.vp.contigs[node] = g.vp.contigs[g.vertex(v)]
        new_g.vp.ab[node] = g.vp.ab[g.vertex(v)]
    for e in g.edges():
        v1 = old_to_new[int(e.source())]
        v2 = old_to_new[int(e.target())]
        e_new = new_g.add_edge(v1, v2)
        new_g.ep.ori[e_new] = g.ep.ori[e]
    new_paths = {}
    for k, p in paths.items():
        new_p = []
        for (v, ori) in p:
            new_p.append((old_to_new[v], ori))
        new_paths[k] = new_p
    return new_g, new_paths

def get_node_abundances(graph, node_dict, aln_file):
    """
    Reads vg alignments and computes node abundance values and edge abundances
    (connection counts per edge).
    Returns node abundance list and connection count dict.
    """
    import json
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


def compress_graph(g, contig_paths, node_abundances):
    """
    Compress the variation graph: merge any non-branching paths into a single
    node and update all node properties accordingly. Returns the updated graph
    and correspondingly updated contig paths.
    """
    # add dummy nodes at sign flips in paths
    dummy_list = []
    for k, path in contig_paths.items():
        prev_ori = path[0][1]
        prev_v = path[0][0]
        for (v, ori) in path:
            if ori != prev_ori:
                # sign flip; add dummy node
                dummy = g.add_vertex()
                dummy_list.append(dummy)
                # e = g.edge(prev_v, v)
                # if e == None:
                #     e = g.edge(v, prev_v)
                #     ori_v = "+" if g.ep.ori[e][0] == "-" else "-"
                #     ori_v_prev = "+" if g.ep.ori[e][1] == "-" else "-"
                #     if e == None:
                #         print("ERROR: missing edge {}->{} for path {}".format(prev_v, v, k))
                # else:
                #     ori_v = g.ep.ori[e][1]
                #     ori_v_prev = g.ep.ori[e][0]
                #
                # e1 = g.add_edge(prev_v, dummy)
                # g.ep.ori[e1] = "{}+".format(g.ep.ori[e][0])
                # e2 = g.add_edge(dummy, v)
                # g.ep.ori[e2] = "+{}".format(g.ep.ori[e][1])
                e1 = g.add_edge(prev_v, dummy)
                g.ep.ori[e1] = "{}+".format(prev_ori)
                e2 = g.add_edge(dummy, v)
                g.ep.ori[e2] = "+{}".format(ori)
            prev_ori = ori
            prev_v = v

    # find simple edges
    simple_edges = []
    simple_adj_in = {}
    simple_adj_out = {}
    for e in g.edges():
        # check subgraph based on vertex orientations as dictated by this edge
        u = e.source()
        v = e.target()
        ori_e = g.ep.ori[e]
        out_degree_u = 0
        in_degree_v = 0
        for e_out in u.out_edges():
            if g.ep.ori[e_out][0] == ori_e[0]:
                out_degree_u += 1
        for e_in in u.in_edges():
            if g.ep.ori[e_in][1] != ori_e[0]:
                out_degree_u += 1
        for e_out in v.out_edges():
            if g.ep.ori[e_out][0] != ori_e[1]:
                in_degree_v += 1
        for e_in in v.in_edges():
            if g.ep.ori[e_in][1] == ori_e[1]:
                in_degree_v += 1

        # if e.source().out_degree() == 1 and e.target().in_degree() == 1:
        if in_degree_v == 1 and out_degree_u == 1:
            assert e.source() != e.target()
            if e.source().out_degree() == 1 and e.target().in_degree() == 1:
                # also take graph structure into account to avoid
                # merging conflicts (e.g. v0 <-(+-)- v1 -(++)-> v2)
                simple_edges.append([e.source(), e.target()])
                simple_adj_out[int(e.source())] = e
                simple_adj_in[int(e.target())] = e

    # build simple paths from simple edges
    def extend_path(p):
        v = int(p[-1])
        if v in simple_adj_out:
            p.append(simple_adj_out[v].target())
            return extend_path(p)
        else:
            return p
    simple_paths = []
    for v, e in simple_adj_out.items():
        if v not in simple_adj_in:
            p = extend_path([e.source(), e.target()])
            simple_paths.append(p)

    # merge simple paths and store new node abundances in graph
    del_list = []
    replacement_dict = {}
    g.vp.ab = g.new_vertex_property('float')
    for v in g.vertices():
        try:
            g.vp.ab[v] = node_abundances[int(v)]
        except KeyError as e:
            if not v in dummy_list:
                print(v)
                print(node_abundances)
            assert v in dummy_list
    for path in simple_paths:
        # print([int(x) for x in path])
        # compute average weighted abundance of path
        total_len = 0
        total_cov = 0
        total_seq = ""
        contigs = []
        for v in path:
            seq = g.vp.seq[v]
            total_seq += seq
            total_len += len(seq)
            total_cov += len(seq) * g.vp.ab[v]
            contigs += g.vp.contigs[v]
        av_cov = total_cov / total_len
        # merge path into start node
        start_node = path[0]
        g.vp.seq[start_node] = total_seq
        g.vp.ab[start_node] = av_cov
        g.vp.contigs[start_node] = list(set(contigs))
        # update outgoing edges
        for v in path[-1].out_neighbors():
            e = g.edge(path[-1], v)
            e_new = g.add_edge(start_node, v)
            # if v == start_node:
            #     print("ERROR: creating self-loop at {}".format(start_node))
            #     print(path)
            #     sys.exit(1)
            g.ep.ori[e_new] = g.ep.ori[e]
        # also add incoming edges:
        # in normal orientation these will be removed as soon as simple
        # nodes are removed, while in opposite orientation they will be kept
        for u in path[-1].in_neighbors():
            if u != start_node:
                e = g.edge(u, path[-1])
                e_new = g.add_edge(u, start_node)
                g.ep.ori[e_new] = g.ep.ori[e]
        del_list += path[1:]
        for v in path[1:]:
            replacement_dict[v] = path[0]
    # add dummy nodes to deletion list
    del_list += dummy_list
    del_list = list(set(del_list))
    # keep track of vertex IDs that were kept, in order to map to new IDs
    keep_list = []
    for v in g.vertices():
        if v not in del_list:
            keep_list.append(int(v))
    # assert keep_list == sorted(keep_list)

    # update contig paths
    old2new_IDs = {x : i  for i, x in enumerate(keep_list)}
    updated_paths = {}
    err_count = 0
    for k, path in contig_paths.items():
        new_path = []
        for (v, ori) in path:
            try:
                v_new = old2new_IDs[int(v)]
                local_ori = ori
            except KeyError:
                try:
                    v_new = old2new_IDs[int(replacement_dict[v])]
                except KeyError:
                    err_count += 1
                    break
                    print("path id", k)
                    print(path)
                    print("v, ori", v, ori)
                    print()
                    print("simple paths:")
                    for p in simple_paths:
                        print([int(x) for x in p])
                    print()
                    print("old2new_IDs:")
                    print(old2new_IDs)
                    sys.exit(1)
                local_ori = ori
            if len(new_path) == 0 or new_path[-1][0] != v_new:
                new_path.append((v_new, ori))
        assert len(new_path) > 0
        updated_paths[k] = new_path
    print("err_count = {}".format(err_count))
    # now remove nodes from simple paths
    # (don't do this earlier because it affects replacement_dict)
    g.remove_vertex(del_list)
    # assert has_simple_paths(g) == False

    return g, updated_paths


def write_gfa(g, gfa_file, paths=None):
    """
    Write graph-tool graph to a gfa format storing nodes, edges, and paths.
    Returns nothing.
    """
    contigs = {} # dict mapping contigs to lists of nodes
    with open(gfa_file, 'w') as f:
        f.write("H\tVN:Z:1.0") # header line GFA 1.0
        for v in g.vertices():
            f.write("\nS\t{}\t{}".format(int(v), g.vp.seq[v]))
            for w in v.out_neighbors():
                e = g.edge(v, w)
                ori = g.ep.ori[e]
                assert len(ori) == 2
                f.write("\nL\t{}\t{}\t{}\t{}\t0M".format(int(v), ori[0],
                                                         int(w), ori[1]))
            for k in g.vp.contigs[v]:
                if k in contigs.keys():
                    contigs[k].append((v, "+"))
                else:
                    contigs[k] = [(v, "+")]

        if paths != None:
            if len(paths) > 0:
                # write paths from path list
                path_info = paths
            else:
                # write paths based on graph traversal
                path_info = contigs
                print("WARNING: using contig paths generated by graph traversal")
            for k, node_list in path_info.items():
                if len(node_list) == 0:
                    print("Skipping length 0 path")
                    continue
                path = ""
                match = ""
                # for v in contigs[k]:
                for (v, ori) in node_list:
                    try:
                        vertex = g.vertex(v)
                        path += str(int(v)) + ori + ','
                        match += str(len(g.vp.seq[vertex]))+'M,'
                    except ValueError:
                        print(k)
                        print(v)
                        print(node_list)
                        print(g.vertex(v))
                        sys.exit(1)
                path = path.rstrip(',')
                match = match.rstrip(',')
                f.write("\nP\t{}\t{}\t{}".format(k, path, match))
    return


def filter_edges(contig_graph, paths, edge_abundances, min_edge_ab, branching_only):
    """
    Filters edges in graph based on a minimal edge abundance,
    then updates graph and contig paths.
    """
    if branching_only:
        edge_del_list = [edge for edge, ab in edge_abundances.items() \
                        if (ab < min_edge_ab and (edge.source().out_degree() > 1
                        or edge.target().in_degree() > 1))]
    else:
        edge_del_list = [edge for edge, ab in edge_abundances.items() \
                        if ab < min_edge_ab]

    contig_blocks = {c : [path] for c, path in paths.items()}
    for e in edge_del_list:
        v1 = e.source()
        v2 = e.target()
        contigs = set(contig_graph.vp.contigs[v1]).intersection(
                    set(contig_graph.vp.contigs[v2]))
        print("Removing edge {} -> {}".format(v1, v2))
        # print("checking contigs {}".format(v1, v2, contigs))
        for c in contigs:
            # check contig path
            split_contig = []
            for path in contig_blocks[c]:
                node, ori = path[0]
                block_start = 0
                next_idx = 1
                for next_node, next_ori in path[1 : ]:
                    if (node, next_node) == (int(v1), int(v2)) \
                            or (next_node, node) == (int(v1), int(v2)):
                        # break contig (TODO: reroute if possible)
                        # TODO: check orientations
                        split_contig.append(path[block_start : next_idx])
                        block_start = next_idx
                    node = next_node
                    next_idx += 1
                split_contig.append(path[block_start : ])
            contig_blocks[c] = split_contig
    # update paths (new identifiers for broken contigs)
    new_paths = {}
    for c, blocks in contig_blocks.items():
        if len(blocks) == 1:
            new_paths[c] = blocks[0]
        else:
            print("Breaking contig {} into {} blocks".format(c, len(blocks)))
            for i, path in enumerate(blocks):
                block_id = "{}_{}".format(c, i)
                new_paths[block_id] = path
                # update contig graph
                for v, ori in path:
                    contig_list = list(contig_graph.vp.contigs[v])
                    contig_list.remove(c)
                    contig_list.append(block_id)
                    contig_graph.vp.contigs[v] = contig_list
    # finally remove edges from graph
    for e in edge_del_list:
        contig_graph.remove_edge(e)
    return contig_graph, new_paths


def read_gfa(gfa_file):
    """
    Reads a graph from a GFA-file and returns graph in gt-format.
    """
    print("Reading graph from {}".format(gfa_file))
    # Define a graph with its vertex properties
    g = Graph(directed=True)
    g.vp.seq = g.new_vertex_property('string')
    g.vp.contigs = g.new_vertex_property('vector<string>')
    g.ep.ori = g.new_edge_property('string')

    # read gfa and add vertices to graph
    node_dict_old2new = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'S':
                # vertex line
                node_id = line[1]
                v = g.add_vertex()
                assert node_id not in node_dict_old2new
                node_dict_old2new[node_id] = int(v)
                seq = line[2].upper()
                if len(seq) == 0:
                    print("WARNING: empty sequence in GFA")
                    print(line)
                g.vp.seq[v] = seq

    # parse through gfa again to add edges and contig paths to graph
    paths = {}
    path_count = 0
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'L':
                # edge line
                v1 = node_dict_old2new[line[1]]
                v2 = node_dict_old2new[line[3]]
                e = g.add_edge(v1, v2)
                g.ep.ori[e] = "{}{}".format(line[2], line[4])
            elif line[0] == 'P':
                # path line
                path_count += 1
                path_name = line[1]
                path_line = line[2].split(',')
                path = [(node_dict_old2new[node_inf[:-1]],
                         node_inf[-1]) for node_inf in path_line]
                if path_name not in paths.keys():
                    paths[path_name] = path
                else:
                    print("WARNING: sequence {} has multiple "
                          "alignments".format(path_name))
                    while path_name in paths.keys():
                        path_name = path_name + "*"
                    paths[path_name] = path
                for (node, ori) in path:
                    g.vp.contigs[node].append(path_name)

    print("Vertex count: {}".format(len(list(g.vertices()))))
    print("Edge count: {}".format(len(list(g.edges()))))
    print("Path count: {}".format(path_count))

    ordered_paths = OrderedDict(sorted(paths.items(), key=lambda x:x[0]))
    return g, ordered_paths, node_dict_old2new


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


def is_cycle(cycle, adj_out):
    """Check if a given cycle still exists according to adjacency list"""
    for i, u in enumerate(cycle):
        if i == len(cycle)-1:
            v = cycle[0]
        else:
            v = cycle[i+1]
        if v not in adj_out[u]:
            return False
    return True


def check_parallel_nodes(g):
    """Returns a list of matching subsequences for parallel nodes."""
    # print("get_dup_nodes")
    dup_list = []
    for v in g.vertices():
        # check out-branches for parallel sequence
        neighbors = list(v.out_neighbors())
        if len(neighbors) <= 1:
            continue
        for i, w1 in enumerate(neighbors):
            for w2 in neighbors[i+1 : ]:
                e1 = g.edge(v, w1)
                e2 = g.edge(v, w2)
                ori1 = g.ep.ori[e1]
                ori2 = g.ep.ori[e2]
                if ori1[0] != ori2[0]:
                    # print("start-orientations don't match")
                    continue
                seq1 = g.vp.seq[w1]
                seq2 = g.vp.seq[w2]
                score = align(seq1, ori1[1], seq2, ori2[1])
                if score > max(len(seq1), len(seq2)) / 2:
                    dup_list += [w1, w2]
        # also check in-branches
        neighbors = list(v.in_neighbors())
        if len(neighbors) <= 1:
            continue
        for i, w1 in enumerate(neighbors):
            for w2 in neighbors[i+1 : ]:
                e1 = g.edge(w1, v)
                e2 = g.edge(w2, v)
                ori1 = g.ep.ori[e1]
                ori2 = g.ep.ori[e2]
                if ori1[1] != ori2[1]:
                    # print("end-orientations don't match")
                    continue
                seq1 = g.vp.seq[w1]
                seq2 = g.vp.seq[w2]
                score = align(seq1, ori1[0], seq2, ori2[0])
                if score > max(len(seq1), len(seq2)) / 2:
                    dup_list += [w1, w2]
    dup_list = list(set(dup_list))
    print("# of nodes with substantial parallel sequence: {}".format(len(dup_list)))
    return dup_list


def align(seq1, ori1, seq2, ori2):
    """Align a given pair of sequences, return alignment blocks."""
    if ori1 == ori2:
        score = pairwise2.align.globalxx(seq1, seq2,
                                         score_only=True)
    elif ori1 == "-":
        score = pairwise2.align.globalxx(rev_comp(seq1), seq2,
                                         score_only=True)
    else:
        score = pairwise2.align.globalxx(seq1, rev_comp(seq2),
                                         score_only=True)
    return score

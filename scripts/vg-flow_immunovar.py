#!/usr/bin/python3

import sys
from numpy import array, zeros, multiply, dot, ceil, where, mod
from gurobipy import *
#sys.path.append('~/tools/vg-flow/scripts/') # use 'export PYTHONPATH=$PYTHONPATH:$vg_flow_dir'
from graph_tool.all import Graph
from graph_tool.topology import is_DAG, all_circuits, topological_sort, shortest_distance, shortest_path
import argparse
# from tqdm import tqdm # progress tracker
import datetime
import time
import math
import subprocess
import copy
from collections import OrderedDict

from graph_functions import read_gfa, get_seq, is_cycle, rev_comp

TIME_LIMIT = 3600
sys.setrecursionlimit(10**6)

__author__ = "Adaptations by Dylan Duchen - Original code by Jasmijn Baaijens"
__license__ = "MIT"

usage = 'Build strains from contigs in variation graph and estimate their abundances.'

def main():
    parser = argparse.ArgumentParser(prog='python opt_edge.py', description=usage)
    parser.add_argument('abundancefile', type=str, help='Node abundance file')
    parser.add_argument('graphfile', type=str, help='GFA file containing the contig variation graph')
    parser.add_argument('--optimization_approach', dest='obj_func', type=int, default=1, help='Optimzation approach (1-7)')
    parser.add_argument('-m', '--min_abundance', dest='min_ab', type=int, required=True, help="Minimal node abundance; nodes below this threshold are considered erroneous and allowed to match any other node.")
    parser.add_argument('-c', '--min_cov', dest='min_cov', type=float, required=True, help='Minimum coverage required per strain')
    parser.add_argument('-d', '--min_depth', dest='min_depth', type=int, default=50, help='Output a list of nodes with sequence depth less than <min_depth>.')
    parser.add_argument('-o', '--out_fasta', dest='fasta', type=str, default='haps.final.fasta', help='Output fasta file with final haplotypes')
    parser.add_argument('-p', '--paths', type=str, help='Read candidate paths from file instead of enumerating them.')
    parser.add_argument('-r', '--reduce_obj', dest='reduce_obj', type=int, default=0, help='Use a reduced objective function by specifying how many times a given combination of paths will be evaluated (reduces runtime and memory usage).')
    parser.add_argument('-t', '--threads', type=int, default=1, help="Set number of threads used for Gurobi.")
    #parser.add_argument('-cb', '--check_branches', action='store_true', help='disable paths that cross branches without contig evidence when available.')
    parser.add_argument('--trim', dest='trim', type=int, default=10, help='number of bases trimmed on either end of contig')
    parser.add_argument('--remove_included_paths', dest='included_paths', type=int, default=1, help='remove overlapping contigs/paths (1=yes, 0=no)')
    parser.add_argument('--max_strains', dest='max_strains', type=int, default=0, help='set an upper bound on the number of strains allowed when using ILP-based reconstruction')
    parser.add_argument('--toboggan', dest='toboggan', type=str, help='toboggan executable')
    parser.add_argument('--ilp', dest='ilp', action='store_true', help='reconstruct strains by solving ILP')
    parser.add_argument('--greedy_mode', dest='mode', default='all')
    parser.add_argument('--careful', dest='careful', action='store_true', help='turn on careful mode to break cycles more extensively. Leads to more fragmented assemblies to avoid misassemblies')
    parser.add_argument('--time_limit', dest='time_limit', default=36000)
    #parser.add_argument('-sp', '--store_paths', action='store_true', help='store all paths and strains in text files.')
    args = parser.parse_args()

    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    print()

    global TIME_LIMIT
    TIME_LIMIT = float(args.time_limit)

    # Flow optimization problem - multiple solutions
    # 1=based on minimizing absolute differences
    # 2= minimizing relative absolute differences - linear approach: abs(abundance - sum_xv) / abundance
    # 3= minimizing relative absolute differences sqrt approach: abs(abundance - sum_xv) / math.sqrt(abundance)
    # 4 or 7= minimize squared error (7 with strain penalty): (abundance - sum_xv) * (abundance - sum_xv)
    # 5 = relative squared error (abundance - sum_xv) * (abundance - sum_xv) / max(0.00001, abundance)

    global obj_func
    obj_func = float(args.obj_func)

    if args.obj_func == 1:
        print(
            "based on minimizing absolute differences"
        )
    if args.obj_func == 2:
        print(
            "minimizing relative absolute differences - linear approach: abs(abundance - sum_xv) / abundance"
        )
    if args.obj_func == 3:
        print(
            "minimizing relative absolute differences sqrt approach: abs(abundance - sum_xv) / math.sqrt(abundance)"
        )
    if args.obj_func == 4:
        print(
            "minimize squared error (7 with strain penalty): (abundance - sum_xv) * (abundance - sum_xv)"
        )
    if args.obj_func == 5:
        print(
            "relative squared error (abundance - sum_xv) * (abundance - sum_xv) / max(0.00001, abundance)"
        )

    if args.ilp and args.max_strains == 0:
        print(
            "ERROR: ILP-mode requires an upper bound on the number of strains. "
            "Use the --max_strains argument to set the maximum number of "
            "strains to be reconstructed."
        )
        sys.exit(1)

    # read node abundances from file
    abundance_list = []
    low_ab_nodes = []
    min_ab = args.min_ab
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
                low_ab_nodes.append(int(node_id))

    # read graph from file
    #graph = load_graph(args.graphfile)
    graph, paths = read_gfa(args.graphfile)[:2]
    nvert = len(list(graph.vertices()))
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

    # trim contigs and get contig adjacency matrix
    graph, paths, adj_out, start_info, vg_to_contigs = build_adj_info(
                graph, paths, nvert, low_ab_nodes, args.trim, args.included_paths)
    write_intermediate_files(graph, paths, adj_out, start_info)
    # check if contig overlap graph is cyclic
    print("Checking graph for cycles...")
    adj_graph, node_to_contig = graphtool_graph(adj_out)
    break_all, break_in, break_out = break_cycles(adj_graph, len(adj_out))
    if len(break_all) == 0:
        print("contig adjacency graph is acyclic")
    else:
        print("breaking cycles at {} points".format(len(break_all)))
        if args.careful:
            break_in = [node_to_contig[v] for v in break_in]
            break_out = [node_to_contig[v] for v in break_out]
        else:
            break_in = [node_to_contig[v] for v in break_all]
            break_out = [node_to_contig[v] for v in break_all]
        new_adj_out = {}
        for v, adj_list in adj_out.items():
            if v in break_out:
                new_adj_out[v] = []
            else:
                new_adj_out[v] = [w for w in adj_list if w not in break_in]
        adj_out = new_adj_out
        adj_graph, node_to_contig = graphtool_graph(adj_out)
    assert is_DAG(adj_graph)

    # num_cycles = 0
    # print("enumerating cycles (if any) in contig graph...")
    # for cycle in all_circuits(adj_graph, unique=True):
    #     num_cycles += 1
    #     cycle_contigs = [node_to_contig[v] for v in cycle]
    #     if not is_cycle(cycle_contigs, adj_out):
    #         continue
    #     broken = False
    #     for i, v in enumerate(cycle_contigs):
    #         if len(adj_out[v]) > 1:
    #             print("removing edges {} -> {}".format(v, adj_out[v]))
    #             for v2 in cycle[i].out_neighbors():
    #                 for v3 in v2.in_neighbors():
    #                     adj_out[node_to_contig[v3]] = []
    #             broken = True
    #             break
    #     if not broken:
    #         v = cycle_contigs[-1]
    #         print("disjoint cycle; removing edge {} -> {}".format(v, adj_out[v]))
    #         adj_out[v] = []
    # if num_cycles == 0:
    #     print("contig adjacency graph is acyclic")
    # else:
    #     print("{} cycles removed from contig adjacency graph".format(num_cycles))

    # cycle_list = cyclic(adj_out)
    # if len(cycle_list) > 0:
    #     # print("contig adjacency graph is cyclic, exiting")
    #     # sys.exit(1)
    #     print("contig adjacency graph contains {} cycles.".format(
    #             len(cycle_list)))
    #     print("breaking cycles")
    #     for cycle in cycle_list:
    #         print("cycle: ", cycle)
    #         broken = False
    #         for v in cycle:
    #             if len(adj_out[v]) > 1:
    #                 print("removing edges {} -> {}".format(v, adj_out[v]))
    #                 adj_out[v] = []
    #                 broken = True
    #                 break
    #         if not broken and is_cycle(cycle, adj_out):
    #             v = cycle[-1]
    #             print("removing edges {} -> {}".format(v, adj_out[v]))
    #             adj_out[v] = []
    # else:
    #     print("contig adjacency graph is acyclic")

    # check branches: if any contigs bridge the branch, ignore all contigs that
    # start/stop at the internal branch node
    discarded = check_branches(graph, adj_out)
    transitive_edges = get_transitive_edges(adj_out)
    for (v1, v2) in set(discarded + transitive_edges):
        adj_out[v1].remove(v2)

    skip_nodes = get_extremity_nodes(graph, args.min_depth)
    print("\nFound {} extremity nodes (depth < {})\n".format(len(skip_nodes),
                                                                args.min_depth))

    # Now solve the minimization problem
    minimization_min_cov = 0
    x, objVal, flow_graph, node_dict = optimize_flow(
            abundance_list, vg_to_contigs, adj_out, skip_nodes,
            args.reduce_obj, minimization_min_cov, args.min_cov,
            args.min_ab, args.threads)
    # write contigs and their subpaths to files for error investigation
    contig_dict = {}
    for i, c in enumerate(adj_out):
        contig_dict[i] = c
    contig_abundances = get_contig_ab(flow_graph, x, contig_dict, node_dict)
    write_intermediate_files(graph, paths, adj_out, start_info, contig_abundances)

    if args.toboggan:
        # Decompose flow into paths using toboggan
        toboggan_graph = write_toboggan_graph(flow_graph, x)
        # Run toboggan
        toboggan_output = "toboggan.out"
        print("\n***Running toboggan to find flow decomposition***\n")
        subprocess.check_call("python {} {} --output {}".format(
                args.toboggan, toboggan_graph, toboggan_output), shell=True)
        # Interpret output
        cg_paths, sol = toboggan_to_paths(toboggan_output, node_dict)
    elif args.ilp:
        # Reconstruct strains using ILP
        print("\n***Reconstruct strains using ILP***")
        cg_paths, sol = optimize_haps(
                args.min_cov, args.max_strains, x, flow_graph, contig_dict,
                args.threads)
    else:
        # Reconstruct strains using greedy algo
        print("\n***Reconstruct strains using greedy approach***")
        if args.mode == "all":
            cg_paths = []
            sol = []
            for mode in ["max_capacity", "min_capacity", "shortest"]:
                graph_copy = copy.deepcopy(flow_graph)
                x_copy = copy.deepcopy(x)
                cg_paths_m, sol_m = greedy_haps(args.min_cov, x_copy, graph_copy,
                                                node_dict, mode)
                for path, value in zip(cg_paths_m, sol_m):
                    if path not in cg_paths:
                        cg_paths.append(path)
                        sol.append(value)
        else:
            cg_paths, sol = greedy_haps(args.min_cov, x, flow_graph,
                                        node_dict, args.mode)

    print("***Processing output***")
    print("{} candidate paths:".format(len(cg_paths)))
    # print(cg_paths)
    # print("contig paths:")
    # print(paths)
    # print(sol)
    print("Contig adjacencies:", adj_out)
    vg_paths, haps, sol = contig_paths_to_nodes(cg_paths, graph, paths,
            start_info, low_ab_nodes, abundance_list, min_ab, sol)
    if not args.ilp:
        # Solve final LP to obtain optimal abundance estimates for greedy paths
        sol, objVal = optimize_abundances(abundance_list, nvert, vg_paths,
            skip_nodes, args.reduce_obj, args.max_strains, args.min_cov,
            args.threads)
    final_abundances = process_output(sol, [], haps, args.min_cov, args.fasta)
    build_genome_graph(graph, vg_paths, final_abundances)
    return

################################################################################


def get_transitive_edges(adj_out):
    # translate out-edges to in-edges
    adj_in = {c : [] for c in adj_out.keys()}
    for c, adj_list in adj_out.items():
        for i in adj_list:
            adj_in[i].append(c)
    # find transitive edges
    transitive_edges = []
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
                    transitive_edges.append((v1, i))
                elif v1 in adj_out[v2]:
                    transitive_edges.append((v2, i))
    return transitive_edges


def get_contig_ab(flow_graph, flow_sol, contig_dict, node_dict):
    abundances = {}
    # print(list(flow_graph.edges()))
    # print(flow_sol)
    assert (len(list(flow_graph.edges())) == len(flow_sol))
    for i, e in enumerate(flow_graph.edges()):
        source_node =  int(e.source())
        if source_node > 1 and source_node%2 == 0:
            contig_id = contig_dict[(source_node-2)/2]
            # print(i, e, flow_sol[i])
            # print(contig_id, source_node, node_dict[contig_id])
            assert node_dict[contig_id][0] == source_node
            abundances[contig_id] = flow_sol[i]
    return abundances


def greedy_haps(min_cov, x, flow_graph, node_dict, mode, skip_overlap_flow=True):
    contig_paths = []
    abundances = []
    edge_index = {}
    flow_graph.ep.weight = flow_graph.new_edge_property('int')
    max_flow = max(x)
    nw = False
    i = 0
    for e in flow_graph.edges():
        edge_index[e] = i
        if mode == 'shortest':
            flow_graph.ep.weight[e] = 1
        elif skip_overlap_flow and int(e.source())%2 == 1:
            # avoid using misleading flow values on overlap edges
            flow_graph.ep.weight[e] = max_flow
            x[i] = max_flow
        elif mode == 'max_capacity':
            flow_graph.ep.weight[e] = -x[i]
            nw = True
        elif mode == 'min_capacity':
            flow_graph.ep.weight[e] = x[i]
        else:
            print("greedy reconstruction mode error, exiting.")
            sys.exit(1)
        i += 1
    inverse_node_dict = {}
    for v, tup in node_dict.items():
        v_in, v_out = tup
        inverse_node_dict[v_in] = v
        inverse_node_dict[v_out] = v
    # extract paths one by one greedily
    source = flow_graph.vertex(0)
    sink = flow_graph.vertex(1)
    # while max(x) >= min_cov:
    assert is_DAG(flow_graph)
    print("mode = {}".format(mode))
    while shortest_path(flow_graph, source, sink)[1]:
        vlist, elist = shortest_path(flow_graph, source, sink,
                        weights=flow_graph.ep.weight,
                        negative_weights=nw)
        freq = min([x[edge_index[e]] for e in elist])
        for e in elist:
            x[edge_index[e]] -= freq
            if x[edge_index[e]] < max(1, min_cov):
                flow_graph.remove_edge(e)
            elif mode == 'max_capacity':
                flow_graph.ep.weight[e] = -x[edge_index[e]]
            elif mode == 'min_capacity':
                flow_graph.ep.weight[e] = x[edge_index[e]]
        if freq >= min_cov:
            contig_path = []
            for idx, v in enumerate(vlist):
                if idx%2 == 1 and v > 1:
                    contig_path.append(inverse_node_dict[v])
            contig_paths.append(contig_path)
            abundances.append(freq)
        # print(shortest_path(flow_graph, source, sink))
    print("{} paths extracted greedily in {} mode".format(len(contig_paths),
        mode))
    return contig_paths, abundances


def solve_contig_ILP(num_strains, min_cov, contig_abundances, graph,
        top_ordering, dist_map, threads):
    nvert = len(list(graph.vertices()))
    max_cov = 1.05*max(contig_abundances)
    m = Model('qp')
    obj = QuadExpr()
    f = m.addVars(list(range(num_strains)), lb=0, ub=max_cov,
        vtype=GRB.CONTINUOUS, name='f')
    V = m.addVars(list(range(num_strains*nvert)), vtype=GRB.BINARY, name='V')
    # add extra variable Z to represent f[s]*V[i,s]
    Z = m.addVars(
        list(range(num_strains*nvert)), lb=0, vtype=GRB.CONTINUOUS, name='Z')

    # add path constraints
    for i, v in enumerate(top_ordering):
        vertex = graph.vertex(v)
        for w in top_ordering[i+1:]:
            if dist_map[vertex][w] > nvert:
                for s in range(0, num_strains):
                    m.addConstr(V[int(v) * num_strains + s] +
                        V[int(w) * num_strains + s] <= 1)
        if vertex.in_degree() > 0:
            for s in range(0, num_strains):
                sum_expr = LinExpr()
                for w in vertex.in_neighbors():
                    sum_expr += V[int(w) * num_strains + s]
                m.addConstr(sum_expr >= V[int(v) * num_strains + s])
        if vertex.out_degree() > 0:
            for s in range(0, num_strains):
                sum_expr = LinExpr()
                for w in vertex.out_neighbors():
                    sum_expr += V[int(w) * num_strains + s]
                m.addConstr(sum_expr >= V[int(v) * num_strains + s])

    # force f=0 for empty strains
    for s in range(0, num_strains):
        max_f = 0
        for v in top_ordering:
            max_f += V[int(v) * num_strains + s]
        m.addConstr(f[s] <= max_f * max(contig_abundances) * 100)

    # # no empty paths
    # for s in range(0, num_strains):
    #     sum_expr = LinExpr()
    #     for v in top_ordering:
    #         sum_expr += V[int(v) * num_strains + s]
    #     m.addConstr(sum_expr >= 1)

    # add dummy variable constraints
    for v in graph.vertices():
        for s in range(0, num_strains):
            m.addConstr(Z[int(v) * num_strains + s] <=
                max_cov * V[int(v) * num_strains + s])
            m.addConstr(Z[int(v) * num_strains + s] <= f[s])
            m.addConstr(Z[int(v) * num_strains + s] >=
                f[s] - max_cov * (1-V[int(v) * num_strains + s]))

    # define the objective
    print('\nDefine the objective function:')
    n_eval = 0
    for v in graph.vertices():
        abundance = contig_abundances[int(v)]
        # sum the calculated abundances of strains
        sum_xv = 0
        for s in range(0, num_strains):
            # sum_xv += f[s] * V[int(v)*num_strains + s]
            sum_xv += Z[int(v) * num_strains + s]
        # compute squared error
        # obj += (abundance - sum_xv) * (abundance - sum_xv)
        obj += (abundance - sum_xv) * (abundance - sum_xv) / max(0.00001, abundance)
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
    m.Params.timeLimit = TIME_LIMIT

    #Minimize the model for the given objective function and constraints
    print("\n*** Phase 1 optimization***\n")
    try:
        m.optimize()
    except:
        return None

    print("Number of solutions = {}".format(m.solcount))
    if m.status == GRB.Status.OPTIMAL or m.status == GRB.Status.TIME_LIMIT:
        f_sol = []
        V_sol = []
        nstrains_test = 0
        for v in m.getVars():
            if 'f' in v.varName:
                f_sol.append(v.x)
            elif 'V' in v.varName:
                V_sol.append(v.x)
        print('\nObjective value: %g' % m.objVal)
        obj_val = m.objVal
        return obj_val, V_sol, f_sol

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
        print('\nNo optimal solution found.')
        return None


def optimize_haps(min_cov, max_strains, flow_sol, flow_graph, node_dict, threads):
    # translate flow graph to contig variation graph
    contig_graph = Graph(directed=True)
    contig_abundances = []
    for v in flow_graph.vertices():
        if int(v) > 1 and int(v)%2 == 0:
            contig_graph.add_vertex()
    idx = 0
    for e in flow_graph.edges():
        source_node = int(e.source())
        target_node = int(e.target())
        if source_node > 1 and source_node%2 == 0:
            contig_abundances.append(flow_sol[idx])
        elif source_node > 1 and target_node > 1:
            new_source = int((source_node-3)/2)
            new_target = int(target_node/2-1)
            contig_graph.add_edge(new_source, new_target)
        idx += 1
    assert len(contig_abundances) == len(node_dict)
    top_ordering = topological_sort(contig_graph)

    # reconstruct haplotypes
    opt_paths = []
    opt_f = []
    opt_obj_val = 0
    #for num_strains in range(1, max_strains+1):
    dist_map = shortest_distance(contig_graph)
    for num_strains in range(max_strains, max_strains+1):
        sol = solve_contig_ILP(num_strains, min_cov, contig_abundances,
                contig_graph, top_ordering, dist_map, threads)
        if sol != None:
            obj_val, V, f = sol
            if obj_val < opt_obj_val or opt_paths == []:
                opt_paths = V
                opt_f = f
                opt_obj_val = obj_val

    # translate solution to paths through contig graph
    cg_paths = []
    opt_abundances = []
    num_strains = len(opt_f)
    nvert = len(contig_abundances)
    for s in range(num_strains):
        path = []
        for i in top_ordering:
            if opt_paths[i*num_strains + s] == 1:
                path.append(node_dict[i])
        if path != []:
            cg_paths.append(path)
            opt_abundances.append(opt_f[s])
    return cg_paths, opt_abundances


def toboggan_to_paths(toboggan_output, node_dict):
    contig_paths = []
    abundances = []
    inverse_node_dict = {}
    for v, tup in node_dict.items():
        v_in, v_out = tup
        inverse_node_dict[v_in] = v
        inverse_node_dict[v_out] = v
    with open(toboggan_output, 'r') as f:
        i = 0
        for line in f:
            if i%2 == 0:
                id_info = line.strip()
            else:
                path = [int(x) for x in line.strip().split(',')]
                contig_path = []
                for idx, x in enumerate(path):
                    if idx%2 == 1 and x > 1:
                        contig_path.append(inverse_node_dict[x])
                        # contig_path.append("p{}".format(int(x/2-1)))
                contig_paths.append(contig_path)
                abundances.append(float(id_info.split()[1]))
            i += 1
    return contig_paths, abundances


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


def contig_paths_to_nodes(cg_paths, g, paths, start_info, low_ab_nodes,
        abundances, min_ab, sol):
    """
    Translates a list of contig paths (through the contig concatenation graph)
    to the corresponding node paths through the contig variation graph.
    Returns a list of node paths and a list of the corresponding haplotype
    sequences.
    """
    print()
    node_paths = []
    haps = []
    for contig_path in cg_paths:
        node_path = []
        idx_v_to_correct = []
        current_idx = 0
        i = 0
        c_prev = contig_path[0]
        if len(contig_path) > 1:
            c = contig_path[1]
            start_index, e_ori = start_info[c_prev][c]
            if e_ori[0] == "+":
                node_path += paths[c_prev]
            else:
                for (c2, ori) in paths[c_prev][::-1]:
                    flipped_ori = "+" if ori == "-" else "-"
                    node_path.append((c2, flipped_ori))
        else:
            node_path += paths[c_prev]
        for c in contig_path[1:]:
            try:
                start_idx, p_ori = start_info[c_prev][c]
            except KeyError:
                print("WARNING: Constraint violation. Skipping path remainder.")
                print(c_prev, c)
                break
            current_idx += (start_idx - i)
            # traverse overlap to find low-abundance nodes
            if p_ori[1] == "+":
                p = paths[c] # contig path [(node, ori)]
            else:
                p = []
                for (c2, ori) in paths[c][::-1]:
                    flipped_ori = "+" if ori == "-" else "-"
                    p.append((c2, flipped_ori))
            assert len(p) > 0
            assert current_idx >= 0
            i = 0
            while current_idx < len(node_path) and i < len(p):
                if node_path[current_idx] != p[i][0]:
                    idx_v_to_correct.append(current_idx)
                current_idx += 1
                i += 1
            # add remaining path sequence
            node_path += p[i:]
            c_prev = c
        # reconstruct the corresponding haplotype sequence
        hap = []
        for v, ori in node_path:
            hap.append(get_seq(g, v, ori))
        # correct low-abundance nodes
        # corrected_path, hap = correct_path(g, idx_v_to_correct, node_path,
        #                 hap, contig_path, abundances, min_ab)
        corrected_path = node_path
        hap = ''.join(hap)
        haps.append(hap)
        node_paths.append(corrected_path)

    # filter out any duplicate paths
    remove_path_idx = []
    for h1 in range(len(haps)):
        hap1 = haps[h1]
        for h2 in range(h1+1, len(haps)):
            hap2 = haps[h2]
            if hap1 in hap2:
                remove_path_idx.append(h1)
                sol[h2] += sol[h1]
            elif hap2 in hap1:
                remove_path_idx.append(h2)
                sol[h1] += sol[h2]
    for idx in sorted(list(set(remove_path_idx)), reverse=True):
        node_paths.pop(idx)
        haps.pop(idx)
        sol.pop(idx)

    print("Reduced #paths from {} to {}.".format(len(cg_paths), len(node_paths)))

    with open("haps.fasta", 'w') as f_haps:
        for i, hap in enumerate(haps):
            f_haps.write(">path{}\n{}\n".format(i, hap))
    return node_paths, haps, sol


def optimize_flow(
        abundance_list, vg_to_contigs, adj_out, skip_nodes, reduce_obj,
        minimization_min_cov, min_cov, min_ab, threads):
    # define new model
    m = Model('lp')
    flow_graph, node_dict, vg_to_fg = build_flow_graph(adj_out, vg_to_contigs)
    # edge variables
    x = m.addVars(list(flow_graph.edges()), lb=0, vtype=GRB.CONTINUOUS, name='x')
    # add additional variables to implement absolute values
    y = m.addVars(vg_to_contigs, lb=0, vtype=GRB.CONTINUOUS, name='y')

    # set constraints for flow conservation
    for v in flow_graph.vertices():
        if int(v) < 2:
            # source/sink node
            continue
        sum_flow = LinExpr()
        for u in v.in_neighbors():
            sum_flow += x[flow_graph.edge(u,v)]
        for w in v.out_neighbors():
            sum_flow -= x[flow_graph.edge(v,w)]
        m.addConstr(sum_flow == 0)

    # define objective
    obj = LinExpr()
    path_combinations_seen = {}
    count_reduce_effect = 0
    for v, edge_list in vg_to_fg.items():
        assert len(edge_list) > 0
        if v in skip_nodes:
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
        sum_xe = 0
        for e in edge_list:
            weight = flow_graph.ep.weight[e]
            sum_xe += weight * x[e]
        abundance = abundance_list[int(v)]
        # if abundance == 0:
        #     continue
        obj += y[v] #abs(abundance - sum_xv)
        # print("node error:", sum_xe-abundance)
        # set constraints on y[v] to obtain absolute value
        m.addConstr(y[v] >= sum_xe - abundance, "y_{}_-".format(v))
        m.addConstr(y[v] >= abundance - sum_xe, "y_{}_+".format(v))

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

    # solve ILP
    print("\n*** Running optimization***\n")
    m.optimize()

    print("Number of solutions = {}".format(m.solcount))
    if m.status == GRB.Status.OPTIMAL:
        # print("\n***Final optimization results:***")
        x_final = []
        for v in m.getVars():
            # print(v)
            if 'x' in v.varName:
                x_final.append(v.x)
        # print('\nObjective value: %g' % m.objVal)
        objVal = m.objVal
        # print("x_final:", x_final)
        return x_final, objVal, flow_graph, node_dict

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


def write_toboggan_graph(flow_graph, x, outfile='toboggan.graph'):
    with open(outfile, 'w') as f:
        f.write("# graph number = 0 name = flow_graph\n")
        f.write("{}\n".format(len(list(flow_graph.vertices()))))
        i = 0
        for e in flow_graph.edges():
            f.write("{} {} {}\n".format(e.source(), e.target(), x[i]))
            i += 1
    return outfile


def build_flow_graph(adj_out, vg_to_contigs):
    # build flow graph
    g = Graph(directed=True)
    vprop = g.new_vertex_property('string')
    g.vp.id = vprop
    eprop = g.new_edge_property('int')
    g.ep.weight = eprop
    # add global source and sink
    source_node = g.add_vertex()
    g.vp.id[source_node] = 's'
    sink_node = g.add_vertex()
    g.vp.id[sink_node] = 't'
    # add 2 nodes for every contig
    node_dict = {}
    for v in adj_out.keys():
        v_in = g.add_vertex()
        v_out = g.add_vertex()
        e = g.add_edge(v_in, v_out) # contig-edge
        g.ep.weight[e] = 1
        node_dict[v] = (v_in, v_out)
    # add edges for overlaps
    for v, adj_list in adj_out.items():
        v_in, v_out = node_dict[v]
        if len(adj_list) == 0:
            e = g.add_edge(v_out, sink_node) # sink-edge
            g.ep.weight[e] = 0
        for w in adj_list:
            w_in, w_out = node_dict[w]
            e = g.add_edge(v_out, w_in) # overlap-edge
            g.ep.weight[e] = -1
    # add edges for source nodes
    for v in adj_out.keys():
        v_in, v_out = node_dict[v]
        if v_in.in_degree() == 0:
            e = g.add_edge(source_node, v_in) # source-edge
            g.ep.weight[e] = 0
    # build map from vg nodes to flow graph edges
    vg_to_edges = {}
    for v, contigs in vg_to_contigs.items():
        assert len(contigs) > 0
        edges = []
        for c in contigs:
            if isinstance(c, str):
                c_in, c_out = node_dict[c]
                e = g.edge(c_in, c_out)
                assert e != None
                edges.append(e)
            else:
                v1 = node_dict[c[0]][1]
                v2 = node_dict[c[1]][0]
                e = g.edge(v1, v2)
                if e != None:
                    edges.append(e)
        vg_to_edges[v] = edges
    return g, node_dict, vg_to_edges


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
            for v, ori in path:
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
            f.write(",".join([str(old_to_new_nodes[v]) + ori for v, ori in p]) + "\t")
            f.write("M,".join([str(node_lengths[old_to_new_nodes[v]]) for v, ori in p]) + "M")
            f.write("\n")
    return


# def read_gfa(gfa_file):
#     """
#     Reads a graph from a GFA-file and returns graph in gt-format.
#     """
#     # Define a graph with its vertex properties
#     g = Graph(directed=True)
#     vprop = g.new_vertex_property('string')
#     g.vp.seq = vprop
#     vprop = g.new_vertex_property('vector<string>')
#     g.vp.contigs = vprop
#
#     # read gfa and add vertices to graph
#     with open(gfa_file, 'r') as f:
#         for line in f:
#             line = line.rstrip('\n').split('\t')
#             if line[0] == 'S':
#                 # vertex
#                 node_id = int(line[1]) - 1
#                 v = g.add_vertex()
#                 seq = line[2].upper()
#                 if len(seq) == 0:
#                     print(line)
#                 g.vp.seq[v] = seq
#
#     # parse through gfa again to add edges and contig paths to graph
#     path_count = 0
#     paths = {}
#     with open(gfa_file, 'r') as f:
#         for line in f:
#             line = line.rstrip('\n').split('\t')
#             if line[0] == 'L':
#                 # edge
#                 assert line[2] == line[4] == "+"
#                 v1 = int(line[1]) - 1
#                 v2 = int(line[3]) - 1
#                 g.add_edge(v1, v2)
#             elif line[0] == 'P':
#                 path_count += 1
#                 contig_id = line[1]
#                 path = line[2].split(',')
#                 unoriented_path = []
#                 old_ori = ''
#                 for node_info in path:
#                     # take care of node orientations
#                     ori = node_info[-1]
#                     if old_ori != '':
#                         # check if orientations change within the path,
#                         # indicating the presence of inversions
#                         assert ori == old_ori
#                     old_ori = ori
#                     node = int(node_info.rstrip(ori)) - 1
#                     unoriented_path.append(node)
#                     # store contig in node
#                     g.vp.contigs[node].append(contig_id)
#                 if ori == "+":
#                     paths[contig_id] = unoriented_path
#                 else:
#                     paths[contig_id] = unoriented_path[::-1]
#     # # save graph in gt format
#     # g.save(graph_name + '.gt')
#     ordered_paths = OrderedDict(sorted(paths.items(), key=lambda x:x[0]))
#     return g, ordered_paths


def write_trimmed_contig_info(graph, paths):
    """Reconstruct contig sequences from their node paths and write to fasta."""
    contig_seqs = {}
    with open("trimmed_contigs.fasta", 'w') as f:
        with open("trimmed_contigs.paths", 'w') as f2:
            for i, (c, path) in enumerate(paths.items()):
                seq = ""
                path_string = ""
                for v, ori in path:
                    seq += get_seq(graph, v, ori)
                    path_string += "{},".format(v)
                f.write(">{}\n{}\n".format(c, seq))
                f2.write(">{}\n{}\n".format(c, path_string[:-1]))
                contig_seqs[c] = seq
    return contig_seqs


def write_intermediate_files(graph, paths, adj_out, start_info, contig_ab=None):
    """Write trimmed contigs, their node paths and the overlap graph to file."""
    contig_seqs = {}
    with open("trimmed_contigs.fasta", 'w') as f:
        with open("trimmed_contigs.paths", 'w') as f2:
            for i, (c, path) in enumerate(paths.items()):
                try:
                    abundance = "{:1f}x".format(contig_ab[c])
                except TypeError as e:
                    abundance = "."
                seq = ""
                path_string = ""
                for v, ori in path:
                    seq += get_seq(graph, v, ori)
                    path_string += "{},".format(v)
                f.write(">{} ab={}\n{}\n".format(c, abundance, seq))
                f2.write(">{}\n{}\n".format(c, path_string[:-1]))
                contig_seqs[c] = seq
    # write contig-adjacency graph to gfa
    with open("trimmed_contigs.gfa", 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        for c, adj_list in adj_out.items():
            f.write("S\t{}\t{}\n".format(c, contig_seqs[c]))
            for w in adj_list:
                start_index, ori = start_info[c][w]
                f.write("L\t{}\t{}\t{}\t{}\t0M\n".format(c, ori[0], w, ori[1]))


def build_adj_info(graph, paths, nvert, low_ab_nodes, trim_size, included_paths):
    """
    Trims contig extremities, removes inclusions due to trimming,
    and finds feasible contig concatenations.
    """
    # trim contigs
    trimmed_list = []
    del_list = []
    for c, path in paths.items():
        if len(path) == 0:
            print("Removing path {} of length 0".format(c))
            del_list.append(c)
            continue
        trim_left = 0
        while len(path) > 0:
            (v, ori) = path[0]
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
            (v, ori) = path[-1]
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
        if trim_left > 0 or trim_right > 0:
            trimmed_list.append(c)
        paths[c] = path

    # remove inclusions (due to trimming)
    # NOTE: very naive implementation, improve by traversing graph once
    print("removing inclusions due to trimming/extension")
    path_strings = {c : ''.join([get_seq(graph, v, ori) for v, ori in path])
                        for c, path in paths.items()}
    # path_strings = {c : "_".join([str(x) for x in path])
    #                     for c, path in paths.items()}
    # print(path_strings)
    # for c1 in trimmed_list:
    #     path1 = path_strings[c1]
    for c1, path1 in path_strings.items():
        if len(path1) == 0:
            del_list.append(c1)
            continue
        for c2, path2 in path_strings.items():
            if c1 == c2:
                continue
            elif path1 == path2 or rev_comp(path1) == path2:
                if c1 < c2:
                    if  included_paths == 1:
                        del_list.append(c1)
                        print("{} == {}".format(c1, c2))
                        print("removing the shorter contig {}".format(c1))
                    elif included_paths == 0:
                        print("{} == {}".format(c1, c2))
                        print("retaining both")
                    else:
                        print("use valid option for --remove_included_paths (0=retain/1=remove)")
            elif path1 in path2 or rev_comp(path1) in path2:
                if included_paths == 1:
                    del_list.append(c1)
                    print("{} in {}".format(c1, c2))
                    print("removing {}".format(c1))
                elif included_paths == 0:
                    print("{} in {}".format(c1, c2))
                    print("Retaining both")
                else:
                    print("use valid option for --remove_included_paths (0=retain/1=remove)")
    del_set = set(del_list)

    # build overlap list per contig (out-edges)
    dup_list = get_dup_nodes(graph)

    # adj_out = {}
    # start_info = {}
    # starts_per_node = {v : [] for v in range(nvert)}
    # for c, path in paths.items():
    #     v = path[0][0]
    #     starts_per_node[v].append(c)
    # for c, path in paths.items():
    #     adj_list, start_indexes = get_adj_out(c, paths, graph, low_ab_nodes,
    #             starts_per_node, dup_list, del_set)
    #     adj_list = [x for x in adj_list if x not in del_set]
    #     adj_out[c] = adj_list
    #     start_info[c] = {x : start_indexes[x] for x in adj_list}
    #     # also add overlap for all vertices contained
    #     for x in adj_list:
    #         start_idx = start_info[c][x]
    #         overlap = path[start_idx:]
    #         for v, ori in overlap:
    #             vg_to_contigs[v].append([c, x])


    # write trimmed contig sequences to file
    write_trimmed_contig_info(graph, paths)
    low_ab_nodes += dup_list
    adj_out, start_info, inc = get_adj_out(paths, graph, low_ab_nodes, del_set)
    # remove any further inclusions found due to low abundance nodes
    del_set.update(inc)
    print("# inclusions found: {}".format(len(del_set)))
    for c in del_set:
        del paths[c]
    # store for every node in the variation graph to which contigs it belongs
    print("update contig info per node")
    vg_to_contigs = {}
    for v in graph.vertices():
        contigs = list(set(graph.vp.contigs[v]) - del_set)
        if contigs:
            vg_to_contigs[v] = contigs

    # keep ordered dict to ensure deterministic behaviour
    print("keep ordered adjacency lists")
    ord_adj_out = OrderedDict(sorted(adj_out.items(), key=lambda x:x[0]))
    print("add overlap info")
    for c, adj_list in adj_out.items():
        for x in adj_list:
            start_idx, ori = start_info[c][x]
            if ori[0] == "+":
                overlap = paths[c][start_idx:]
            else:
                overlap = paths[c][::-1][start_idx:]
            for v, ori in overlap:
                vg_to_contigs[v].append([c, x])
    return graph, paths, ord_adj_out, start_info, vg_to_contigs


def check_branches(graph, adj_out):
    """
    Check branches in overlap graph: if any branching edge is supported by a
    clique of size 3, discard all overlaps that are not supported by such a
    clique. Returns a list of contig overlaps to be discarded during path
    enumeration.
    """
    print("check_branches")

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
    print("# of duplicate nodes: {}".format(len(dup_list)))
    return dup_list


def get_adj_out(paths, g, low_ab_nodes, del_set):
    """
    Returns an adjacency list & info of all feasible contig concatenations.
    """
    adj_out = {c : [] for c in paths.keys() if c not in del_set}
    adj_in = {c : [] for c in paths.keys() if c not in del_set}
    start_info = {c : {} for c in adj_out.keys()}
    minimap2_overlaps = "overlaps.minimap2.paf"
    input = "trimmed_contigs.fasta"
    subprocess.check_call(
        "minimap2 -X -x ava-ont --end-bonus=100 {0} {0} > {1}".format(
            input, minimap2_overlaps),
        shell=True
    )

    inclusions = []
    with open(minimap2_overlaps, 'r') as f:
        print("processing overlaps")
        i = 0
        for line in f:
            i += 1
            if i%1000 == 1:
                print("overlaps status: line", i)
            line = line.strip('\n').split('\t')[:12]
            pos1 = int(line[2])
            ori = line[4]
            pos2 = int(line[7])
            if ori == "+":
                if pos1 > pos2:
                    c1 = line[0]
                    c2 = line[5]
                else:
                    c2 = line[0]
                    c1 = line[5]
            else:
                c1 = line[0]
                c2 = line[5]
            if c1 == c2 or c1 in del_set or c2 in del_set:
                continue
            start_index, edge_ori = check_overlap(paths[c1], paths[c2],
                                                      ori, low_ab_nodes)
            if start_index < 0:
                # dismiss overlap (alignment not good enough)
                continue
            elif c2 in adj_out[c1]:
                # overlap already exists
                continue
            # elif len(paths[c1]) - start_index >= len(paths[c2]):
            #     # c2 included in c1
            #     # print("inclusion {} in {}".format(c2, c1))
            #     inclusions.append(c2)
            #     continue
            # elif start_index == 0:
            #     # c1 included in c2
            #     # print("inclusion {} in {}".format(c1, c2))
            #     inclusions.append(c1)
            #     continue
            adj_out[c1].append(c2)
            adj_in[c2].append(c1)
            start_info[c1][c2] = (start_index, edge_ori)
    print("initial adj_out ready: {} edges".format(sum([len(x) for c, x in adj_out.items()])))

    def traverse_from_node(c):
        # print("traverse from node {}".format(c))
        # print("#contig orientations:", len(contig_orientations))
        ori_found = False
        # check incoming edges for vertex orientations
        for c_in in adj_in[c]:
            if ori_found:
                break
            try:
                ori_c_in = contig_orientations[c_in]
            except KeyError as e:
                continue
            # determine contig orientation
            start_index, edge_ori = start_info[c_in][c]
            if ori_c_in == edge_ori[0]:
                contig_orientations[c] = edge_ori[1]
            else:
                ori_c = "+" if edge_ori[1] == "-" else "-"
                contig_orientations[c] = ori_c
            ori_found = True
        # check outgoing edges for vertex orientations
        for c_out in adj_out[c]:
            if ori_found:
                break
            try:
                ori_c_out = contig_orientations[c_out]
                ori_found = True
            except KeyError as e:
                continue
            # determine contig orientation
            start_index, edge_ori = start_info[c][c_out]
            if ori_c_out == edge_ori[1]:
                contig_orientations[c] = edge_ori[0]
            else:
                ori_c = "+" if edge_ori[0] == "-" else "-"
                contig_orientations[c] = ori_c
        if not ori_found:
            contig_orientations[c] = "+"
        neighbors = set(adj_in[c] + adj_out[c])
        for c2 in neighbors:
            if c2 not in contig_orientations:
                traverse_from_node(c2)
        return
    # assign contig orientations
    print("assign contig orientations")
    contig_orientations = {}
    for c in paths.keys():
        if c in del_set:
            continue
        elif c not in contig_orientations:
            traverse_from_node(c)
    # update edge orientations & start info
    updated_adj_out = {c : [] for c in paths.keys() if ((c not in del_set) and (c not in inclusions))}
    updated_start_info = {c : {} for c in updated_adj_out.keys()}
    conflict_list = []
    print("update edge orientations in adj_out")
    num_keys = len(paths.keys())
    for i, c in enumerate(paths.keys()):
        # print("{} / {}".format(i, num_keys))
        if c in del_set or c in inclusions:
            continue
        neighbors = adj_out[c]
        ori_c = contig_orientations[c]
        # print("contig {} - {} neighbors".format(c, len(neighbors)))
        for c2 in neighbors:
            # print("neighbor {}".format(c2))
            if c2 in inclusions:
                continue
            ori_c2 = contig_orientations[c2]
            try:
                start_index, edge_ori = start_info[c][c2]
            except KeyError as e:
                print("start_info not found?!")
                print(c, c2)
                print(start_info[c])
                sys.exit(1)
            # print("check edge")
            if edge_ori[0] == ori_c and edge_ori[1] == ori_c2:
                # print("edge ok")
                updated_adj_out[c].append(c2)
                updated_start_info[c][c2] = (start_index, edge_ori)
                continue
            elif edge_ori[0] != ori_c and edge_ori[1] != ori_c2:
                # flip edge
                # print("flipping edge")
                new_ori = "{}{}".format(ori_c2, ori_c)
                new_index = len(paths[c2]) - len(paths[c]) + start_index
                updated_adj_out[c2].append(c)
                updated_start_info[c2][c] = (new_index, new_ori)
            else:
                # conflicting edge; disconnect node outwards
                # print("conflicting edge")
                conflict_list.append([c, c2])
                continue
    print("final adj_out ready: {} edges".format(sum([len(x) for c, x in updated_adj_out.items()])), flush=True)
    print("# conflicting edges resolved: {}".format(len(conflict_list)))
    return updated_adj_out, updated_start_info, inclusions


def check_overlap(path1, path2, p_ori, low_ab_nodes):
    """
    Check overlaps for a given pair of paths given a suggestion orientation;
    return the start index in the first path if overlap, otherwise return -1.
    """
    if p_ori == "-":
        flipped_path2 = []
        for (v, ori) in path2[::-1]:
            flipped_ori = "+" if ori == "-" else "-"
            flipped_path2.append((v, flipped_ori))
        edge_ori = "+-"
        start_index = align_paths(path1, flipped_path2, low_ab_nodes)
        if start_index == -1:
            flipped_path1 = []
            for (v, ori) in path1[::-1]:
                flipped_ori = "+" if ori == "-" else "-"
                flipped_path1.append((v, flipped_ori))
            edge_ori = "-+"
            start_index = align_paths(flipped_path1, path2, low_ab_nodes)
    else:
        edge_ori = "++"
        start_index = align_paths(path1, path2, low_ab_nodes)
    return start_index, edge_ori


def align_paths(path1, path2, low_ab_nodes):
    """
    Align a given pair of paths; return the start index in the
    first path if overlap, otherwise return -1.
    """
    start_index = -1
    # find the first high abundance node in path2 to search for in path1
    for i, (v, ori) in enumerate(path2):
        if v not in low_ab_nodes:
            search_node = v
            search_index = i
            break
        return start_index # no overlap
    # search for path overlap
    for i, (u, ori_u) in enumerate(path1):
        if u == search_node:
            start_index = i - search_index # take low ab nodes into account
            j = search_index
            overlap = path1[i : ]
            assert len(overlap) > 0
            for (v, ori_v) in overlap:
                if j >= len(path2):
                    # inclusion
                    break
                elif (v, ori_v) != path2[j] and v not in low_ab_nodes \
                                          and path2[j][0] not in low_ab_nodes:
                    # not an acceptable overlap
                    start_index = -1
                    break
                j += 1
            if start_index >= 0:
                # overlap found
                return start_index
    return start_index


# def get_adj_out(c, paths, g, low_ab_nodes, starts_per_node, dup_nodes, del_set):
#     """
#     Returns an adjacency list of all feasible concatenations for a given contig.
#     """
#     path = paths[c]
#     # print(path)
#     nvert = len(list(g.vertices()))
#     start_node = path[0][0]
#     active_contigs = set(g.vp.contigs[g.vertex(start_node)])
#     rescued_contigs = set()
#     start_indexes = {}
#     assert len(active_contigs) > 0
#     contig_seq = [get_seq(g, v, ori) for v, ori in path]
#     # find out-neighbors
#     for i, (node, ori) in enumerate(path):
#         old_active_contigs = active_contigs
#         # update active contig set
#         if node in low_ab_nodes + dup_nodes:
#             assert len(path) > 1
#             if i > 0:
#                 prev_node = path[i-1][0]
#                 alt_nodes = list(g.vertex(prev_node).out_neighbors())
#             else:
#                 next_node = path[i+1][0]
#                 alt_nodes = list(g.vertex(next_node).in_neighbors())
#             for alt in alt_nodes:
#                 active_contigs |= set(starts_per_node[alt])
#         else:
#             current_contigs = set(g.vp.contigs[g.vertex(node)])
#             starts = set(starts_per_node[node])
#             active_contigs = active_contigs.intersection(current_contigs)
#             # print("current contigs: ", current_contigs)
#             # print("starts: ", starts)
#             # print("active contigs: ", active_contigs)
#             active_contigs |= starts
#         # for any discarded active contigs, check if sequence is truly different
#         for c_old in old_active_contigs - active_contigs:
#             if c_old in del_set:
#                 continue
#             idx = start_indexes[c_old]
#             overlap = ''.join(contig_seq[idx:])
#             assert len(overlap) > 0
#             overlap_idx = 0
#             keep_c = True
#             overlap2 = ""
#             for u, ori_u in paths[c_old]:
#                 if not keep_c:
#                     break
#                 u_seq = get_seq(g, u, ori_u)
#                 overlap2 += u_seq
#                 for nuc in u_seq:
#                     if overlap_idx >= len(overlap):
#                         break
#                     elif nuc == 'N' or overlap[overlap_idx] == 'N':
#                         overlap_idx += 1
#                         continue
#                     elif nuc != overlap[overlap_idx]:
#                         keep_c = False
#                         break
#                     overlap_idx += 1
#             if keep_c:
#                 rescued_contigs.add(c_old)
#                 # print("rescued contig overlap {} -> {}".format(c, c_old))
#         # keep track of start points
#         for x in active_contigs:
#             if x not in start_indexes:
#                 start_indexes[x] = i
#     assert c in active_contigs
#     active_contigs |= rescued_contigs
#     active_contigs.remove(c)
#     # TODO: remove paths that have insertions in their overlap
#
#     return list(active_contigs), start_indexes


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


def break_cycles(graph, num_contigs):
    """
    Traverse graph through DFS and returns a list of node IDs corresponding to
    cycle nodes to be broken.
    """
    break_all = []
    break_in = []
    break_out = []
    visited = [False] * num_contigs
    def cycle_dfs(node, path):
        local_path = path.copy()
        visited[int(node)] = True
        local_path[int(node)] = True
        for neighbor in node.out_neighbors():
            if visited[int(neighbor)] == False:
                cycle_dfs(neighbor, local_path)
            elif local_path[int(neighbor)]:
                break_all.append(node)
                for v in neighbor.in_neighbors():
                    break_out.append(int(v))
                for v in neighbor.out_neighbors():
                    break_in.append(int(v))
        return

    for v in graph.vertices():
        if v.in_degree() == 0:
            path = [False] * num_contigs
            cycle_dfs(v, path)

    for v in graph.vertices():
        if not visited[int(v)]:
            path = [False] * num_contigs
            cycle_dfs(v, path)
    assert sum(visited) == num_contigs
    break_all = list(set(break_all))
    break_in = list(set(break_in))
    break_out = list(set(break_in))
    return break_all, break_in, break_out


def graphtool_graph(adj_out):
    """Return a list of cycles given the adjacency lists."""
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
            if c2 not in adj_out.keys():
                print("ERROR: {} not in adj_out while neighbor of {}".format(c2, c))
                print(adj_out)
            graph.add_edge(contig_to_node[c], contig_to_node[c2])
    # if is_DAG(graph):
    #     return [], {}
    # else:
    #     print("enumerating all cycles in contig graph...")
    #     cycle_iterator = all_circuits(graph, unique=True)
    #     # print("converting cycle nodes to contig IDs and store in list")
    #     # cycles = []
    #     # for cycle in cycle_iterator:
    #     #     cycles.append([node_to_contig[v] for v in cycle])
    return graph, node_to_contig


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
                f_out.write(">path{} {}x freq={:.3f}\n".format(l, cov, cov/total_cov))
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


def optimize_abundances(a, nvert, paths, skip_nodes, reduce_obj, max_strains,
                        min_cov_final, threads):
    """
    Defines Gurobi minimization problem and then applies the LP solver.
    Returns the optimal haplotype abundances and the corresponding objective
    value.
    """
    t1_start = time.perf_counter()
    t2_start = time.process_time()
    #obj_func = 1
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
            m.addConstr(x_binary[i] >= (x[i])/max_x)
        # define total strain count
        sum_x_binary = LinExpr()
        for i in range(npaths):
            sum_x_binary += x_binary[i]
        # bound the number of strains
        m.addConstr(sum_x_binary <= max_strains)

    # Store paths in P: p_ij = 1 if node i contains path j
    # print('\nSave for every node which paths are passing through:')
    for i in range(npaths):
        for v, ori in paths[i]:
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
    # print('\nDefine the objective function:')
    n_eval = 0
    for v in range(nvert):
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
    m.Params.timeLimit = TIME_LIMIT
    #m.Params.InfUnbdInfo = 1
    #m.Params.PreQLinearize = 1 # default = automatic

    #Minimize the model for the given objective function and constraints
    print("\n*** Phase 1 optimization***\n")
    m.optimize()

    print("Number of solutions = {}".format(m.solcount))
    if m.status == GRB.Status.OPTIMAL or m.status == GRB.Status.TIME_LIMIT:
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
        t1_stop = time.perf_counter()
        t2_stop = time.process_time()
        print("\nAbundance estimation completed")
        print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
        print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
        return(x_final, objVal)

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


if __name__ == '__main__':
    sys.exit(main())
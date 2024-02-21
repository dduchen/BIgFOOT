#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import json
from tqdm import tqdm
from graph_tool.all import *

def main():
    parser = argparse.ArgumentParser(prog='check_node_abundance.py',
                description='Compare node abundances to ground truth.')
    parser.add_argument('-g', '--graph', dest='graph', type=str,
                required=True, help='Graph in .gt format')
    parser.add_argument('-x', '--extremities', dest='extremity_file', type=str,
                help='Text file containing graph extremity nodes')
    parser.add_argument('-t', '--total_cov', dest='total_cov', type=float,
                required=True, help='Total coverage (average) for the data set')
    parser.add_argument('-a', '--ab_file', dest='ab_file', type=str,
                required=True, help='Input abundance file')
    parser.add_argument('-s', '--strain_file', dest='strain_file', type=str,
                required=True, help='Ground truth strains (.fastq)')
    parser.add_argument('-f', '--freq_file', dest='freq_file', type=str,
                required=True, help='Ground truth strain frequencies')
    parser.add_argument('-o', '--out_file', dest='out_file', type=str,
                default='true_abundance.txt', help='Write true abundances here')
    parser.add_argument('-vg', dest='vg', type=str, default="vg",
                help="Path to vg executable")
    args = parser.parse_args()

    graph_name = args.graph.rstrip('.gt')
    graph = load_graph(args.graph)
    aln_name = graph_name + ".strains_aln"

    # make sure that ground truth strains are in fastq format
    if not args.strain_file[-1] == 'q':
        print("\nERROR: ground truth strains need to be in fastq format.\n")
        sys.exit(1)

    # align true strains to graph
    if os.path.exists(aln_name + '.filtered.json'):
        aln_file = aln_name + '.filtered.json'
        print("\nAlignment file found ---> skipping strain alignment step.")
    else:
        print()
        aln_file = vg_map_strains(graph_name, args.strain_file, aln_name, args.vg)

    # determine strains per node
    strains_per_node, true_paths, missing_count, nodes = get_strains_map(graph, aln_file)

    # determine true node abundances
    strain_frequencies = read_freq_file(args.freq_file)
    print(strain_frequencies)
    true_node_abundances = {}
    with open(args.out_file, 'w') as f:
        for node, strains in strains_per_node.items():
            total_freq = sum([strain_frequencies[s] for s in strains])
            node_cov = total_freq * args.total_cov
            #node_id = str(int(node)-1)
            node_id = str(node)
            true_node_abundances[node_id] = node_cov
            f.write("{}:{:.1f}\n".format(node_id, node_cov))

    # read in the computed node abundance values
    computed_node_abundances = {}
    with open(args.ab_file, 'r') as f:
        renamed_node = 0
        for line in f:
            [node_id, abundance] = line.rstrip().split(':')
            computed_node_abundances[str(renamed_node)] = float(abundance)
            renamed_node += 1

    # compare abundances
    false_node_count, false_nodes = compare_node_abundances(computed_node_abundances,
            true_node_abundances, args.total_cov, args.extremity_file, nodes)

    print("\n\n\n")
    for node in false_nodes:
        print(node)
        v  = graph.vertex(node)
        contigs = graph.vp.strain[v]
        #print(contigs)
    print("\n\n\n")
    #print("False node count: {}".format(false_node_count))
    print("Missing node count: {} (including extremities)".format(missing_count))
    print()
    return


def compare_node_abundances(computed_node_abundances, true_node_abundances,
        total_cov, extremity_file, nodes):
    skip_nodes = []
    if extremity_file:
        with open(extremity_file, 'r') as f:
            for line in f:
                node = line.rstrip()
                skip_nodes.append(node)
    error_list = []
    rel_error_list = []
    false_node_coverage = []
    underestimated = []
    overestimated = []
    false_node_count = 0
    false_nodes = []
    error_file = open("abundance_errors.txt", 'w')
    for node_id, abundance in computed_node_abundances.items():
        if node_id in skip_nodes:
            continue
        error = true_node_abundances[node_id] - abundance
        if true_node_abundances[node_id] > 0:
            rel_error = abs(error)/true_node_abundances[node_id]*100
        else:
            seq = nodes[int(node_id)]
            print("false:", node_id, abundance, seq)
            false_node_count += 1
            false_node_coverage.append(error)
            false_nodes.append(node_id)
            continue
        error_list.append(error)
        rel_error_list.append(rel_error)
        error_file.write("{}:\t{:.1f}\t{:.1f}\t{:.1f}\n".format(
            node_id, error, float(abundance), true_node_abundances[node_id]))
        if error < 0:
            overestimated.append(true_node_abundances[node_id])
        elif error > 0:
            underestimated.append(true_node_abundances[node_id])
    error_file.close()

    abs_error_list = [abs(x) for x in error_list]
    average_error = sum(abs_error_list)/len(abs_error_list)
    rel_average_error = average_error/total_cov*100
    average_rel_error = sum(rel_error_list)/len(rel_error_list)
    print("\naverage error: {:.1f} ({:.1f}% of total coverage)".format(average_error, rel_average_error))
    print("average relative error: {:.1f}%".format(average_rel_error))
    print("max relative error: {:.1f}%".format(max(rel_error_list)))
    print("smallest error: {:.1f}".format(min(abs_error_list)))

    print("\n#underestimated: {}".format(len(underestimated)))
    print("#overestimated: {}".format(len(overestimated)))
    print("average underestimated abundance: {:.0f}".format(
                                    sum(underestimated)/len(underestimated)))
    print("average overestimated abundance: {:.0f}".format(
                                    sum(overestimated)/len(overestimated)))

    error_bound_10perc = 0.1*total_cov
    error_bound_1perc = 0.01*total_cov
    error_bound_01perc = 0.001*total_cov
    count_10perc = sum([1 if x<=error_bound_10perc else 0 for x in abs_error_list])
    count_1perc = sum([1 if x<=error_bound_1perc else 0 for x in abs_error_list])
    count_01perc = sum([1 if x<=error_bound_01perc else 0 for x in abs_error_list])

    print("\nFalse node count: {} out of {} evaluated".format(false_node_count, len(error_list)))
    if false_node_count > 0:
        print("average false node coverage: {:.1f}x\n".format(sum(false_node_coverage)/len(false_node_coverage)))
    # print("#nodes within 10 perc error range: {}".format(count_10perc))
    # print("#nodes within 1 perc error range: {}".format(count_1perc))
    # print("#nodes within 0.1 perc error range: {}".format(count_01perc))

    error_bins = [-10.0, -8.0, -6.0, -4.0, -2.0, -1.0, -0.0, 0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    count_per_bin = [0 for x in error_bins]
    perfect_count = 0
    for error in error_list:
        if error > 0:
            for i in range(1, 1+int(len(error_bins)/2)):
                if error >= error_bins[-i]/100*total_cov:
                    break
            count_per_bin[-i] += 1
        elif error < 0:
            for i in range(int(len(error_bins)/2)):
                if error <= error_bins[i]/100*total_cov:
                    break
            count_per_bin[i] += 1
        else:
            perfect_count += 1
    print("#perfect estimates: {}".format(perfect_count))
    for i in range(len(error_bins)):
        if i <= len(error_bins)/2-1:
            print("#errors <= {} percent coverage: {}".format(error_bins[i], count_per_bin[i]))
        else:
            print("#errors >=  {} percent coverage: {}".format(error_bins[i], count_per_bin[i]))
    print()
    return false_node_count, false_nodes

def vg_map_strains(graph_name, strains_fastq, aln_name, vg, build_index=True):
    kmer_size = 11
    if build_index:
        threads = 1
        subprocess.check_call("mkdir -p tmp", shell=True)
        # convert graph to gfa
        subprocess.check_call(
            "~/git_repos/vg_haps/scripts/gt_gfa_converter.py {0}.gt {0}.gfa".format(graph_name),
            shell=True
        )
        # build graph in vg format
        print("building graph in vg format...")
        subprocess.check_call( vg +
            " view -F -v {0}.gfa > {0}.vg".format(graph_name), shell=True)

        # map reads (global alignment)
        print("indexing graph...")
        #input_check = input("Is your current $TMPDIR large enough to store intermediate gcsa files?")
        subprocess.check_call( "export TMPDIR=tmp && " + vg +
            " index -x {0}.xg -g {0}.gcsa -k {1} -X 2 -t {2} {0}.vg".format(
            graph_name, kmer_size, threads),
            shell=True)
    else:
        print("assuming graph is already indexed!")

    print("mapping reads to graph...")
    subprocess.check_call( vg +
        " map -f {0} -x {1}.xg -g {1}.gcsa -k {2} > {3}.gam".format(
        strains_fastq, graph_name, kmer_size, aln_name),
        shell=True)

    # filter for primary alignments
    print("filter alignments...")
    subprocess.check_call( vg +
        " filter -r 0.90 -s 2 -fu {0}.gam > {0}.filtered.gam".format(aln_name),
        shell=True)

    # output json
    print("convert gam to json format")
    subprocess.check_call( vg +
        " view -a {0}.filtered.gam > {0}.filtered.json".format(aln_name),
        shell=True)
    return "{0}.filtered.json".format(aln_name)


def get_strains_map(graph, aln_file):
    # read the graph
    nodes = {}
    gfa_id_to_real = {}
    gfa_id = 1
    for node in graph.vertices():
        node_id = node
        seq = graph.vp.seq[node]
        nodes[node_id] = seq
        gfa_id_to_real[gfa_id] = node_id
        gfa_id += 1
    # read the alignments and map node IDs to strains
    strains_per_node = {node : [] for node in nodes}
    true_paths = {}
    missing_node_count = 0
    print("\nProcessing alignments...\n")
    with open(aln_file, 'r') as aln_json:
        for line in aln_json:
            aln = json.loads(line)
            seq_id = aln["name"].split()[0]
            seq = aln["sequence"]
            aln_len = 0
            aln_len2 = 0
            aln_path = []
            diff_count = 0
            insertions = []
            deletions = []
            reconstructed_seq = ""
            try:
                path = aln["path"]
                mapping = path["mapping"]
            except KeyError:
                # sequence unmapped
                print("strain {} was not aligned (!)".format(seq_id))
                continue
            #offset = 0
            oldrank = 0
            for node_info in mapping:
                position = node_info["position"]
                node = position["node_id"]
                node_id = gfa_id_to_real[int(node)]
                try:
                    is_reverse = position["is_reverse"]
                    is_reverse = True
                    node_seq = revcomp(nodes[node_id])
                except KeyError:
                    is_reverse = False
                    node_seq = nodes[node_id]
                node_len = len(node_seq)
                rank = int(node_info["rank"])
                assert rank == oldrank + 1
                oldrank = rank
                edit = node_info["edit"]
                for aln_piece in edit:
                    try:
                        diff_length = len(aln_piece["sequence"])
                    except KeyError:
                        diff_length = 0
                    diff_count += diff_length
                    try:
                        from_len = int(aln_piece["from_length"])
                    except KeyError:
                        from_len = 0
                    try:
                        to_len = int(aln_piece["to_length"])
                    except KeyError:
                        to_len = 0
                    if from_len > to_len:
                        deletions.append(from_len - to_len)
                        missing_node_count += 1
                        #print("from_len > to_len for node {}".format(node_id))
                    elif from_len < to_len:
                        insertions.append(-from_len + to_len)
                        missing_node_count += 1
                        #print("from_len < to_len for node {}".format(node_id))
                    else:
                        missing_node_count += int(diff_length > 0)
                        if len(insertions) == 0 or insertions[-1] != '-':
                            insertions.append('-')
                        if len(deletions) == 0 or deletions[-1] != '-':
                            deletions.append('-')
                    aln_len2 += min(from_len, to_len)
                try:
                    offset = int(position["offset"])
                except KeyError:
                    offset = 0
                aln_len += node_len - offset
                if reconstructed_seq == "":
                    reconstructed_seq += node_seq[offset:]
                else:
                    reconstructed_seq += node_seq[:len(node_seq)-offset]
                strains_per_node[node_id].append(seq_id)
                aln_path.append(node_id)
            true_paths[seq_id] = aln_path
            print("Strain {} has length {} and is aligned over {} positions with {} differences (ins: {}, del: {})".format(seq_id, len(seq), aln_len2, diff_count, insertions, deletions))
            print(len(seq), len(reconstructed_seq))
            if insertions[0] != '-':
                seq = seq[insertions[0]:]
            if insertions[-1] != '-':
                seq = seq[:-insertions[-1]]
            print(seq == reconstructed_seq)
    #print("Missing node count: {}".format(missing_node_count))
    return strains_per_node, true_paths, missing_node_count, nodes


def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    revcomp = "".join(complement.get(base, base) for base in reversed(seq))
    assert len(seq) == len(revcomp)
    return revcomp


def read_freq_file(freq_file):
    strain_frequencies = {}
    with open(freq_file, 'r') as f:
        for line in f:
            [strain_id, freq] = line.rstrip('\n').split(':')
            strain_frequencies[strain_id] = float(freq)
    return strain_frequencies


if __name__ == '__main__':
    sys.exit(main())

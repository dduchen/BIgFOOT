#!/usr/bin/python3

import sys, os
import argparse
import subprocess
import time

__author__ = "Jasmijn Baaijens"
__license__ = "MIT"

usage = 'Sort contigs for variation graph reconstruction.'

def main():
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('contigfile', type=str, help='Fasta file containing input contigs')
    parser.add_argument('outfile', type=str, help='Sorted fasta file')
    parser.add_argument('-t', dest='threads', type=int, default=8, help='Maximum number of allowed threads to be used in parallel')
    parser.add_argument('-e', dest='sfo_err', type=float, default=0, help='Maximum relative hamming distance (%) for SFO')
    parser.add_argument('-p', dest='pident', type=float, default=95, help='Maximum relative edit distance (%) for minimap2')
    parser.add_argument('-m', dest='min_len', type=int, default=30, help='Minimal contig overlap length')
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    global v_print
    if args.verbose:
        v_print = print
    else:
        v_print = lambda *a: None  # do-nothing function

    global_t1_start = time.perf_counter()
    global_t2_start = time.process_time()

    overlaps = "overlaps"
    sort_contigs(
        args.contigfile,
        args.outfile,
        overlaps,
        args.threads,
        args.sfo_err,
        args.pident,
        args.min_len,
        args.quick,
        args.verbose
    )
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    print("sort_contigs.py completed")
    print("Total elapsed time: {:.1f} seconds".format(t1_stop-global_t1_start))
    print("Total CPU process time: {:.1f} seconds".format(t2_stop-global_t2_start))
    print()
    return

def sort_contigs(input, output, overlaps, threads, sfo_err,
                 pident, min_len, quick_mode, verbose):
    # 0. remove old sorted_contigs.fasta and its index files
    if os.path.isfile(output):
        subprocess.check_call("rm {}*".format(output), shell=True)

    # 1. find approximate overlaps between input contigs
    t1_start = time.perf_counter()
    t2_start = time.process_time()
    v_print("Searching for approximate overlaps between contigs...")
    find_overlaps(input, overlaps, threads, sfo_err,
                  pident, min_len, quick_mode)
    # report runtime
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    v_print("Overlaps ready!")
    v_print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
    v_print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
    v_print()

    # 2. store overlap lengths in adjacency matrix (overlap graph)
    t1_start = time.perf_counter()
    t2_start = time.process_time()
    v_print("Building adjacency matrix...")
    contig_list, contig_ids, contig_idx_dict = read_fasta(input)
    ncontigs = len(contig_list)
    adj_mat = [[0 for i in range(ncontigs)] for j in range(ncontigs)]
    sfo_overlaps = "{}.sfo.tsv".format(overlaps)
    minimap_overlaps = "{}.minimap2.rf.tsv".format(overlaps)
    adj_mat = update_adj_mat(sfo_overlaps, adj_mat, contig_idx_dict)
    adj_mat = update_adj_mat(minimap_overlaps, adj_mat, contig_idx_dict)
    # report runtime
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    v_print("Adjacency matrix ready!")
    v_print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
    v_print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
    v_print()

    # 3. iteratively add contigs by decreasing overlap lengths
    t1_start = time.perf_counter()
    t2_start = time.process_time()
    if quick_mode:
        v_print("Sorting contigs in quick mode")
        ordering = quick_ordering(contig_list, adj_mat)
    else:
        v_print("Sorting contigs in regular mode")
        ordering = deduce_ordering(contig_list, adj_mat)
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    v_print("Sorting completed")
    v_print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
    v_print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
    v_print()

    # 4. store contigs in new order, convert sequences to uppercase (if necessary)
    t1_start = time.perf_counter()
    t2_start = time.process_time()
    new_fasta = output
    with open(new_fasta, 'w') as f:
        for idx in ordering:
            id = contig_ids[idx]
            seq = contig_list[idx]
            f.write('>{}\n{}\n'.format(id, seq))
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    v_print("Sorted contigs written to {}".format(new_fasta))
    v_print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
    v_print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
    v_print()
    return


def find_overlaps(input, output, threads, sfo_err, pident, min_len, quick_mode):
    # rust-overlaps: all overlaps within 5% hamming distance
    if not quick_mode:
        try:
            subprocess.check_call(
                "rust-overlaps -i -r -w {} {} {}.sfo.tsv {} {}".format(
                threads, input, output, sfo_err, min_len),
                shell=True
            )
        except subprocess.CalledProcessError as e:
            print(e.output)
            print("Skipping rust-overlaps, relying solely on minimap2 for sorting.")
            # create empty overlaps file to prevent downstream errors
            subprocess.check_call("> {}.sfo.tsv".format(output), shell=True)
    else:
        subprocess.check_call("> {}.sfo.tsv".format(output), shell=True)
    # minimap2: heuristic search for overlaps using minimizers)
    subprocess.check_call(
        "minimap2 -X -x ava-ont --end-bonus=100 {0} {0} > {1}.minimap2.paf".format(input, output),
        shell=True
    )
    PAF2SFO(output, min_len, pident)
    return


def PAF2SFO(output, min_len, pident):
    infile = "{}.minimap2.paf".format(output)
    outfile = "{}.minimap2.rf.tsv".format(output)
    too_short_count = 0
    too_div_count = 0
    overlap_count = 0
    with open(outfile, 'w') as f1:
        with open(infile, 'r') as f2:
            for line in f2:
                [qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length, qual] = line.strip('\n').split('\t')[:12]
                if int(length) < min_len:
                    too_short_count += 1
                    continue
                if int(matchcount)/float(length) < pident/100.0:
                    too_div_count += 1
                    continue
                # reformat overlap info
                idA = qseqid
                idB = sseqid
                ori = 'N' if qori == '+' else 'I'
                if ori == 'N':
                    OHA = int(qstart) - int(sstart)
                    OHB = int(slen) - int(sstart) - (int(qlen) - int(qstart))
                else:
                    OHA = int(qstart) - (int(slen) - int(send))
                    OHB = int(send) - (int(qlen) - int(qstart))
                if OHA >= 0:
                    OLA = min(int(qlen) - OHA, int(slen))
                else:
                    OLA = min(int(slen) + OHA, int(qlen))
                OLB = OLA
#                if int(idA) > int(idB):
                if idA > idB:
                    # swap order such that id1 < id2
                    idA, idB = idB, idA
                    if ori == 'N':
                        OHA *= -1
                        OHB *= -1
                    else:
                        # swap orientations such that id1 sequence is forward
                        OHA, OHB = OHB, OHA
                mismatch = int(length) - int(matchcount)
                assert mismatch >= 0
                sfo_line = '\t'.join([idA, idB, ori, str(OHA), str(OHB), str(OLA), str(OLB), str(mismatch)]) + '\n'
                f1.write(sfo_line)
                overlap_count += 1
    print("minimap overlaps shorter than {}bp: {}".format(min_len, too_short_count))
    print("minimap overlaps with less than {} percent identity: {}".format(pident, too_div_count))
    print("total minimap overlaps found: {}".format(overlap_count))
    return


def update_adj_mat(overlaps, adj_mat, contig_idx_dict):
    total_count = 0
    with open(overlaps, 'r') as f:
        for overlap in f:
            total_count += 1
            overlap = overlap.strip()
            # store maximal overlap length in adjacency matrix
            [idA, idB, ori, OHA, OHB, OLA, OLB, mismatch] = overlap.split('\t')
            if idA == idB:
                # ignore self-overlaps
                continue
            len = max(int(OLA), int(OLB))
            idxA = contig_idx_dict[idA]
            idxB = contig_idx_dict[idB]
            if adj_mat[int(idxA)][int(idxB)] < len:
                adj_mat[int(idxA)][int(idxB)] = len
                adj_mat[int(idxB)][int(idxA)] = len
    print("Total number of overlaps: {}".format(total_count))
    return adj_mat


def read_fasta(fasta):
    contig_list = []
    contig_ids = []
    contig_idx_dict = {}
    idx = 0
    seq = ""
    with open(fasta, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line[0] == '>':
                id = line.split()[0].lstrip('>')
                contig_idx_dict[id] = idx
                contig_ids.append(id)
                if idx > 0:
                    contig_list.append(seq)
                idx += 1
                seq = ""
            else:
                seq += line
        contig_list.append(seq)
    return contig_list, contig_ids, contig_idx_dict


def quick_ordering(contig_list, adj_mat):
    """Sort contigs by longest overlaps --- without sorting neighbors"""
    ncontigs = len(contig_list)
    ordering = []
    remaining = [1 for i in range(ncontigs)]
    while sum(remaining) > 0:
        # find longest remaining contig
        max_len = 0
        longest = -1
        for i, r in enumerate(remaining):
            if r == 1 and len(contig_list[i]) > max_len:
                max_len = len(contig_list[i])
                longest = i
        assert longest >= 0
        # ignore overlaps between added elements, add longest remaining contig
        for i in ordering:
            adj_mat[i][longest] = 0
            adj_mat[longest][i] = 0
        ordering.append(longest)
        remaining[longest] = 0
        # add neighbors
        for i, r in enumerate(remaining):
            if r == 1 and adj_mat[i][longest] > 0:
                for j in ordering:
                    adj_mat[j][i] = 0
                    adj_mat[i][j] = 0
                ordering.append(i)
                remaining[i] = 0
    return ordering


def deduce_ordering(contig_list, adj_mat):
    """Sort contigs by longest overlaps"""
    ncontigs = len(contig_list)
    ordering = []
    remaining = [1 for i in range(ncontigs)]
    while sum(remaining) > 0:
        # find longest remaining contig
        max_len = 0
        longest = -1
        for i, r in enumerate(remaining):
            if r == 1 and len(contig_list[i]) > max_len:
                max_len = len(contig_list[i])
                longest = i
        assert longest >= 0
        # ignore overlaps between added elements, add longest remaining contig
        for i in ordering:
            adj_mat[i][longest] = 0
            adj_mat[longest][i] = 0
        ordering.append(longest)
        remaining[longest] = 0
        # iteratively add neighbors
        while sum([(1-remaining[i]) * sum(adj_mat[i]) for i in range(ncontigs)]) > 0:
            # select neighbor c with longest overlap
            max_val = 0
            max_idx = -1
            for row_idx, row in enumerate(adj_mat):
                if remaining[row_idx]:
                    continue
                for col_idx, col in enumerate(row):
                    if col > max_val:
                        max_val = col
                        max_idx = col_idx
            assert max_idx >= 0
            assert max_idx not in ordering
            c = max_idx
            # update ordering, remaining, and adjacencies
            for j in ordering:
                adj_mat[j][c] = 0
                adj_mat[c][j] = 0
            ordering.append(c)
            remaining[c] = 0
    return ordering


if __name__ == '__main__':
    sys.exit(main())

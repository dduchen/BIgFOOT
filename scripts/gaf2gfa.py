#!/usr/bin/env python3

import sys, os
import argparse

# GAF:
# Col	Type	Description
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start (0-based; closed)
# 4	int	Query end (0-based; open)
# 5	char	Strand relative to the path: "+" or "-"
# 6	string	Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
# 7	int	Path length
# 8	int	Start position on the path (0-based)
# 9	int	End position on the path (0-based)
# 10	int	Number of residue matches
# 11	int	Alignment block length
# 12	int	Mapping quality (0-255; 255 for missing)


# GFA paths:
# Column	Field	Type	Regexp	Description
# 1	RecordType	Character	P	Record type
# 2	PathName	String	[!-)+-<>-~][!-~]*	Path name
# 3	SegmentNames	String	[!-)+-<>-~][!-~]*	A comma-separated list of segment names and orientations
# 4	Overlaps	String	\*|([0-9]+[MIDNSHPX=])+	Optional comma-separated list of CIGAR strings


def main():
    parser = argparse.ArgumentParser(description="Add paths from GAF to GFA.")
    parser.add_argument('gaf', type=str, help="GAF alignment file")
    parser.add_argument('gfa', type=str, help="GFA output file to which paths are appended")
    args = parser.parse_args()

    write_gfa_from_gaf(args.gaf, args.gfa)
    return


def write_gfa_from_gaf(gaf_file, gfa_file):
    """Read graph in GAF format and convert to GFA format"""
    gfa = open(gfa_file, 'w')
    with open(gaf_file, 'r') as gaf:
        for line in gaf:
            line = line.rstrip().split()
            name = line[0]
            # get path: a comma-separated list of segment names and orientations
            id = ''
            ori = ''
            path = []
            for char in line[5]:
                if char in '<>':
                    if ori != '':
                        assert id != ''
                        path.append('{}{}'.format(id, ori))
                    id = ''
                    if char == '<':
                        ori = '-'
                    else:
                        ori = '+'
                else:
                    id += char
            path.append('{}{}'.format(id, ori))
            segments = ','.join(path)
            # TODO - get overlaps: a comma-separated list of CIGAR strings
            overlaps = ''
            gfa_line = ["P", name, segments, overlaps]
            gfa.write('\t'.join(gfa_line) + '\n')
    gfa.close()
    return


if __name__ == "__main__":
    sys.exit(main())

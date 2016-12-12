#!/usr/bin/env python
#
# "Reference implementation" of a non-parallel Smith-Waterman algorithm
#
# Created on 09/03/2016
# by Arpad Goretity
#

import sys
import operator
import itertools
import re
import argparse
import numpy


def smith_waterman(seq1, seq2, gap_penalty, score_fn):
    resultMat = numpy.zeros((len(seq1) + 1, len(seq2) + 1), dtype=numpy.int32)

    for i in range(len(seq1) + 1):
        for j in range(len(seq2) + 1):
            if i == 0 and j == 0:
                resultMat[i, j] = 0
            elif i == 0:
                resultMat[i, j] = 0
            elif j == 0:
                resultMat[i, j] = 0
            else:
                resultMat[i, j] = numpy.max([
                    resultMat[i, j - 1] + gap_penalty,
                    resultMat[i - 1, j] + gap_penalty,
                    resultMat[i - 1, j - 1] + score_fn(seq1[i - 1], seq2[j - 1]),
                    0
                ])

    return numpy.array(map(lambda row: row[1:], resultMat[1:]))

def read_n_seqs(f):
    return int(f.readline())

def read_sizes(f):
    return map(int, filter(len, re.split('\\s+', f.readline())))

def read_sequences(f):
    nums = map(int, filter(len, re.split('\\s+', f.read())))
    return list(itertools.izip_longest(*[iter(nums)] * 2))

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Generate and dump sequences and dump scores for the parallelized Smith-Waterman algorithm')
    argparser.add_argument('--scoring-offset', type=int, required=True)
    argparser.add_argument('--gap-penalty', type=int, required=True)
    args = argparser.parse_args(sys.argv[1:])

    n_seqs    = read_n_seqs(sys.stdin)
    seq_sizes = read_sizes(sys.stdin)
    seqs      = read_sequences(sys.stdin)

    def dihedral_score(d1, d2):
        abs_dphi = abs(d1[0] - d2[0])
        abs_dpsi = abs(d1[1] - d2[1])
        dphi = min(abs_dphi, 65536 - abs_dphi)
        dpsi = min(abs_dpsi, 65536 - abs_dpsi)
        return args.scoring_offset - (dphi * dphi + dpsi * dpsi)

    ver_start = 0

    for i in range(0, n_seqs - 1):
        ver_size  = seq_sizes[i]
        ver_end   = ver_start + ver_size
        ver_seq   = seqs[ver_start:ver_end]
        hor_start = ver_end

        sys.stdout.write('#{}.\t'.format(i))

        for j in range(i + 1, n_seqs):
            hor_size = seq_sizes[j]
            hor_end  = hor_start + hor_size
            hor_seq  = seqs[hor_start:hor_end]

            H = smith_waterman(ver_seq, hor_seq, args.gap_penalty, dihedral_score)
            sys.stdout.write(' {}'.format(numpy.max(H)))

            hor_start = hor_end

        ver_start = ver_end
        sys.stdout.write('\n')

    sys.stdout.write('\n')

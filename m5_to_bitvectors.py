import sys
import os
import argparse
import jitu
from dotmap import DotMap

def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference_file", type=str, required=True)
    parser.add_argument("-t", "--transcript", type=str, required=True)
    parser.add_argument("-i", "--input_file", type=str, required=True)
    parser.add_argument("-o", "--output_file", type=str, required=True)
    parser.add_argument("--edge_tolerance", type=int, default=25,
                        help="Allow NA only within this many nts from each end")
    parser.add_argument("--max_ones", type=int, default=None,
                        help="Optional hard cap on number of 1s per read")
    parser.add_argument("--max_ones_frac", type=float, default=None,
                        help="Optional hard cap on fraction of 1s per covered positions")
    args = parser.parse_args()
    return args

def get_state_vector(Ref, d):
    bitD = {i: '.' for i in range(len(Ref))}

    beg = int(d.tStart)
    BLOCK = list(zip(d.tAlignedSeq, d.matchPattern, d.qAlignedSeq))

    if d.tStrand == '-':
        BLOCK = BLOCK[::-1]

    for t, tick, q in BLOCK:
        if t != '-':
            if t == q:
                bitD[beg] = '0'
            else:
                bitD[beg] = '1'
            beg += 1

    return [bitD[i] for i in range(len(Ref))]

def passes_na_filter(state_list, edge_tolerance=25):
    n = len(state_list)
    for i, x in enumerate(state_list):
        if x == '.':
            if i < edge_tolerance:
                continue
            if i >= n - edge_tolerance:
                continue
            return False
    return True

def passes_one_filter(state_list, max_ones=None, max_ones_frac=None):
    ones = sum(1 for x in state_list if x == '1')
    covered = sum(1 for x in state_list if x in ('0', '1'))

    if max_ones is not None and ones > max_ones:
        return False

    if max_ones_frac is not None and covered > 0:
        if ones / covered > max_ones_frac:
            return False

    return True

if __name__ == "__main__":
    args = handler()

    if not os.path.exists(args.input_file):
        print(f"Input file not found: {args.input_file}")
        exit()

    if not os.path.exists(args.reference_file):
        print(f"Reference file not found: {args.reference_file}")
        exit()

    ref_path = args.reference_file
    m5file = args.input_file
    out_bitFile = args.output_file

    tube, seqD = jitu.getTubeD(ref_path)
    assert args.transcript in tube, "check the transcript name in reference file: " + args.transcript

    head = [
        'qName', 'qLength', 'qStart', 'qEnd', 'qStrand',
        'tName', 'tLength', 'tStart', 'tEnd', 'tStrand',
        'score', 'numMatch', 'numMismatch', 'numIns', 'numDel',
        'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq'
    ]

    kept = 0
    dropped_na = 0
    dropped_ones = 0

    with open(m5file) as inp, open(out_bitFile, 'w') as outf:
        for line in inp:
            A = line.strip().split()
            d = DotMap(dict(zip(head, A)))

            if args.transcript != d.tName:
                continue

            Ref = seqD[d.tName]
            state_list = get_state_vector(Ref, d)

            if not passes_na_filter(state_list, edge_tolerance=args.edge_tolerance):
                dropped_na += 1
                continue

            if not passes_one_filter(
                state_list,
                max_ones=args.max_ones,
                max_ones_frac=args.max_ones_frac
            ):
                dropped_ones += 1
                continue

            state_na = ['NA' if x == '.' else x for x in state_list]
            outf.write('\t'.join([d.qName] + state_na) + '\n')
            kept += 1

    print(f"Kept reads: {kept}")
    print(f"Dropped for internal NA: {dropped_na}")
    print(f"Dropped for too many 1s: {dropped_ones}")
    print("DONE")

#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import tempfile
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import jitu

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
fold2dotbracket_path = os.path.join(SCRIPT_DIR, "fold2dotbracketFasta.py")
rnaConvert_path = os.path.join(SCRIPT_DIR, "rnaConvert.py")


def generate_constraints(ref, profile):
    state = []
    for pos, (base, bit) in enumerate(zip(ref, profile), 1):
        if bit == '1':
            state.append([str(pos), base, '0'])
        else:
            state.append([str(pos), base, '-1'])
    return state


def worker_process(index, profile, profile_count, reference_sequence,
                   fold2dotbracket_path, rnaConvert_path):
    """
    Runs contrafold + conversion steps for a single profile.
    Returns a tuple (index, dotbracket, post_content, element_line, profile, profile_count).
    If an error occurs, raises a RuntimeError (will surface to the main process).
    """
    tmpdir = tempfile.mkdtemp(prefix=f"bit_{index}_")
    try:
        state = generate_constraints(reference_sequence, profile)
        constraint_temp = os.path.join(tmpdir, "constraints.bpseq")
        with open(constraint_temp, "w") as outf:
            for i, base, tick in state:
                outf.write(f"{i}\t{base}\t{tick}\n")

        fold_out = os.path.join(tmpdir, "out.fold")
        post_out = os.path.join(tmpdir, "out.post")

        p = subprocess.run(
            ["contrafold", "predict", "--constraints", constraint_temp,
             "--parens", fold_out, "--posteriors", "0.0", post_out],
            capture_output=True, text=True
        )
        if p.returncode != 0:
            raise RuntimeError(f"Contrafold failed for index {index}: {p.stderr.strip()}")

        db_out = os.path.join(tmpdir, "out.db")
        p2 = subprocess.run(
            ["python", fold2dotbracket_path, "--input_file", fold_out,
             "--tag", f"bit_{index}", "--output_file", db_out],
            capture_output=True, text=True
        )
        if p2.returncode != 0:
            raise RuntimeError(f"Fold2dotbracket failed for index {index}: {p2.stderr.strip()}")

        with open(db_out) as f:
            dotbracket = f.read().strip()

        txt_filename = os.path.join(tmpdir, f"bit_{index}.txt")
        p3 = subprocess.run(
            ["python", rnaConvert_path, db_out, "-T", "element_string",
             "--force", "--to-file", "--filename", txt_filename],
            capture_output=True, text=True
        )
        if p3.returncode != 0:
            raise RuntimeError(f"rnaConvert failed for index {index}: {p3.stderr.strip()}")

        element_file = txt_filename + "001.element_string"
        if not os.path.exists(element_file):
            raise RuntimeError(f"Expected element file not found for index {index}: {element_file}")

        element_line = None
        with open(element_file) as inp:
            for line_num, line in enumerate(inp, 1):
                if line_num == 3:
                    element_line = line.strip()
                    break
        if element_line is None:
            raise RuntimeError(f"Element line not found (3rd line) for index {index}")

        with open(post_out) as h:
            post_content = h.read().strip()

        return (index, dotbracket, post_content, element_line, profile, profile_count)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bit_file", type=str, required=False,
                        default="REF.bit")
    parser.add_argument("-r", "--reference_file", type=str, default="REF.fasta")
    parser.add_argument("-t", "--transcript", type=str, default="REFERENCE")
    parser.add_argument("-s", "--size_file", type=str, default=None,
                        help="Optional explicit size filename (if omitted will use bit_basename.size)")
    parser.add_argument("-w", "--workers", type=int, default=None,
                        help="Number of worker processes (defaults to number of CPUs)")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = handler()

    tube, seqD = jitu.getTubeD(fastaFiler=args.reference_file)
    reference_sequence = seqD[args.transcript]

    unique_bits = {}
    profileD = defaultdict(int)
    tot = 0
    with open(args.bit_file) as inp:
        for line in inp:
            A = line.strip().split('\t')
            if A and len(A) > 1:
                bits = ''.join(['1' if b == '1' else '.' for b in A[1:]])
                unique_bits[bits] = bits.count('1')
                profileD[bits] += 1
                tot += 1

    print()
    print(f"Unique bit vectors: {len(unique_bits)} from total: {tot}")

    ordered_profiles = sorted(unique_bits.items(), key=lambda x: x[1], reverse=True)

    bit_basename = os.path.splitext(os.path.basename(args.bit_file))[0]
    merged_db_path = f"{bit_basename}.db"
    merged_txt_path = f"{bit_basename}.txt"
    merged_el_path = f"{bit_basename}.element_string"
    size_filename = args.size_file if args.size_file else f"{bit_basename}.size"

    n_jobs = len(ordered_profiles)
    max_workers = args.workers if args.workers else (os.cpu_count() or 2)
    print(f"There are now {max_workers} lab members getting {n_jobs} unique structures folded.")

    results_by_index = {}
    futures = []
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        for i, (profile, num_mutations) in enumerate(ordered_profiles, start=1):
            pass

    seen = defaultdict(int)
    filtered_job_list = []
    for i, (profile, num_mutations) in enumerate(ordered_profiles, start=1):
        assert len(reference_sequence) == len(profile), "Reference length mismatch!"
        state = generate_constraints(reference_sequence, profile)
        ticks_string = ''.join([x[-1] for x in state])
        seen[ticks_string] += 1
        if seen[ticks_string] > 1:
            continue
        filtered_job_list.append((i, profile, profileD[profile]))

    results_by_index = {}
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        future_to_index = {}
        for (i, profile, pcount) in filtered_job_list:
            fut = exe.submit(worker_process, i, profile, pcount, reference_sequence,
                             fold2dotbracket_path, rnaConvert_path)
            future_to_index[fut] = i

        for fut in tqdm(as_completed(future_to_index), total=len(future_to_index),
                        desc="Time until we get taco bell"):
            try:
                res = fut.result()
                idx = res[0]
                results_by_index[idx] = res
            except Exception as e:
                print(f"\nERROR in worker: {e}", file=sys.stderr)
                raise

    with open(merged_db_path, "w") as merged_db, \
         open(merged_txt_path, "w") as merged_txt, \
         open(merged_el_path, "w") as merged_el, \
         open(size_filename, "w") as outS:

        for (i, profile, pcount) in filtered_job_list:
            if i not in results_by_index:
                raise RuntimeError(f"No result for index {i} — worker likely failed.")
            _, dotbracket, post_content, element_line, profile_str, profile_count = results_by_index[i]

            bit_prefix = f"bit_{i}"
            merged_db.write(dotbracket + "\n")
            merged_txt.write(f"{bit_prefix}\t{post_content}\n")
            merged_el.write(f"{bit_prefix}\t{element_line}\n")
            outS.write(f"{bit_prefix}\t{profile_count}\t{profile_str}\n")

    print("DONE")

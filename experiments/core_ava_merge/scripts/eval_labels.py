"""Evaluate an arm against the tr_substring record labels (base_c8804b1 seed records).

usage: eval_labels.py ARM [ARM ...]

For each labelled seed record (tr_substring/out/record_classes.tsv) the arm's
record carrying the child's allele is found: the record with the same
CHROM/POS/SVTYPE/SVLEN if it carries one of the base record's HG002
consensuses, else the nearest record (<= 2 kb, same type) carrying one of them.
Reported per class: matched, joined to the transmitting parent's consensus
(right), joined only to the other parent's consensus (wrong, when Q100 says
that parent does not transmit), consistent in the arm (strict scorer).
"""

import csv
import gzip
import os
import subprocess
import sys
from bisect import bisect_left
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/home/mayv_c/development/svp_mendel")
PY = "/home/mayv_c/development/svirlpool/.claude/worktrees/strip-error-rate/.pixi/envs/dev/bin/python"
LABELS = ROOT / "tr_substring/out/record_classes.tsv"
SAMPLES = ("HG002", "HG003", "HG004")


def annotate(arm: str) -> Path:
    out = ROOT / "proto" / f"{arm}.mndl.vcf"
    if out.exists():
        return out
    subprocess.run(
        [PY, "-m", "svirlpool.analysis.mendelian_consistency", "-i", str(ROOT / "arms" / arm / "family.vcf.gz"),
         "-o", str(ROOT / "proto" / f"{arm}.mndl.txt"), "-p", str(ROOT / "inputs/giab.ped"), "--vcf", str(out),
         "--min-size", "50", "--passonly"],
        check=True, capture_output=True,
        env={**os.environ, "PATH": str(Path(PY).parent) + ":" + os.environ.get("PATH", "")},
    )
    return out


def read_vcf(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    recs = []
    samples = None
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            t = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
            fmt = t[8].split(":")
            gts = {}
            mndl = None
            for s, v in zip(samples, t[9:], strict=True):
                d = dict(zip(fmt, v.split(":"), strict=False))
                gts[s] = d.get("GT", ".")
                if s == "HG002":
                    mndl = d.get("MNDL")
            recs.append({
                "chrom": t[0], "pos": int(t[1]), "id": t[2], "filter": t[6],
                "svtype": info.get("SVTYPE"), "svlen": abs(int(info.get("SVLEN", "0"))),
                "cons": info.get("CONSENSUSIDs", "").split(","), "gts": gts, "mndl": mndl,
            })
    return recs


def alt(gt: str) -> bool:
    return any(a not in ("0", ".") for a in gt.replace("|", "/").split("/"))


def main():
    labels = list(csv.DictReader(open(LABELS), delimiter="\t"))
    base = read_vcf(ROOT / "tr_substring/data/family.mndl.vcf.gz")
    base_by_key = {(r["chrom"], r["pos"], r["svtype"], r["svlen"]): r for r in base}
    for arm in sys.argv[1:]:
        recs = read_vcf(annotate(arm))
        by_key = defaultdict(list)
        by_cons = defaultdict(list)
        for r in recs:
            by_key[(r["chrom"], r["pos"], r["svtype"], r["svlen"])].append(r)
            for c in r["cons"]:
                by_cons[c].append(r)
        tally = defaultdict(Counter)
        for lab in labels:
            key = (lab["chrom"], int(lab["pos"]), lab["svtype"], int(lab["svlen"]))
            b = base_by_key.get(key)
            k = (lab["cls"], lab["klass"])
            tally[k]["n"] += 1
            if b is None:
                tally[k]["base_missing"] += 1
                continue
            child_cons = [c for c in b["cons"] if c.startswith("HG002:")]
            cand = [r for r in by_key.get(key, []) if set(r["cons"]) & set(child_cons)]
            if not cand:
                near = [
                    r for c in child_cons for r in by_cons.get(c, [])
                    if r["svtype"] == lab["svtype"] and r["chrom"] == lab["chrom"] and abs(r["pos"] - key[1]) <= 2000
                ]
                cand = sorted(near, key=lambda r: abs(r["pos"] - key[1]))
            if not cand:
                tally[k]["lost"] += 1
                continue
            r = cand[0]
            tally[k]["matched"] += 1
            if r["filter"] == "PASS" and r["mndl"] == "consistent":
                tally[k]["consistent"] += 1
            tp = lab["t_parent"]
            other = {"HG003": "HG004", "HG004": "HG003"}.get(tp)
            samples_in = {c.split(":")[0] for c in r["cons"]}
            base_samples = {c.split(":")[0] for c in b["cons"]}
            if tp and tp in samples_in:
                tally[k]["joined_right"] += 1
                if tp not in base_samples:
                    tally[k]["joined_right_new"] += 1
                if alt(r["gts"].get(tp, ".")):
                    tally[k]["right_parent_alt"] += 1
            haps = set(lab["haps"].split(",")) if lab["haps"] else set()
            transmitting = {"MAT": "HG004", "PAT": "HG003"}
            wrong_parent = None
            if haps and len(haps) == 1:
                hp = transmitting[next(iter(haps))]
                wrong_parent = "HG003" if hp == "HG004" else "HG004"
            if wrong_parent and wrong_parent in samples_in and wrong_parent not in base_samples:
                tally[k]["joined_wrong_new"] += 1
                hp = "HG004" if wrong_parent == "HG003" else "HG003"
                if hp not in samples_in:
                    tally[k]["wrong_only_new"] += 1
            if other and other in samples_in and tp not in samples_in:
                tally[k]["joined_other_only"] += 1
        print(f"== {arm}")
        cols = ["n", "matched", "lost", "consistent", "joined_right", "joined_right_new", "right_parent_alt", "joined_other_only", "joined_wrong_new", "wrong_only_new"]
        print("cls\tklass\t" + "\t".join(cols))
        for k in sorted(tally):
            t = tally[k]
            print(f"{k[0]}\t{k[1]}\t" + "\t".join(str(t[c]) for c in cols))


if __name__ == "__main__":
    main()

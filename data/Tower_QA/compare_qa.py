#!/usr/bin/env python3
"""Compare two TowerQA2 outputs and check whether one is a subset of the other.

Works with <outTag>_final_badtowers.txt (key = towerID) and <outTag>_final_badruns.txt (key = runId):
the key is the first column, lines starting with '#' are comments.
Also reads bad-tower lists made of '# <category> ::' headers followed by comma-separated tower IDs
(e.g. the Isobar P20ic EMCal QA list); a tower then counts once, whatever category it is in.

Usage
  python3 compare_qa.py A_final_badtowers.txt B_final_badtowers.txt     # two files of the same kind
  python3 compare_qa.py A_final_badruns.txt   B_final_badruns.txt
  python3 compare_qa.py TowerQA2_quant TowerQA2_sigma                    # two outTags: compares both bad towers and bad runs
  python3 compare_qa.py TowerQA2_final_badtowers.txt P20ic_badtowers.txt # with a category-list file
"""
import os
import re
import sys


def read(path):
    """Return {key: info} for every entry.

    Two formats are understood:
      - TowerQA2 output: one entry per line, key = first column, info = rest of the line ('#' = comment)
      - a list file with '# <category> ::' header lines followed by comma-separated tower IDs
        (e.g. the Isobar P20ic EMCal QA bad towers list); info = the categories the tower is listed in
    """
    with open(path) as f:
        lines = f.read().splitlines()
    if any(l.startswith("#") and "::" in l for l in lines):          # category-list format
        entries, cat = {}, None
        for l in lines:
            if l.startswith("#"):
                cat = l[1:].split("::")[0].strip() if "::" in l else None
                continue
            if cat is None:
                continue
            for tok in re.findall(r"\d+", l):
                entries[tok] = (entries[tok] + "; " if tok in entries else "") + cat
        return entries
    entries = {}
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split(None, 1)
        entries[cols[0]] = cols[1] if len(cols) > 1 else ""
    return entries


def sort_key(k):
    return (0, int(k)) if k.isdigit() else (1, k)


def compare(path_a, path_b, what):
    a, b = read(path_a), read(path_b)
    sa, sb = set(a), set(b)
    only_a, only_b, both = sorted(sa - sb, key=sort_key), sorted(sb - sa, key=sort_key), sa & sb

    name_a, name_b = os.path.basename(path_a), os.path.basename(path_b)
    print(f"=== {what} ===")
    print(f"  A = {name_a}: {len(sa)}")
    print(f"  B = {name_b}: {len(sb)}")
    print(f"  in both: {len(both)},  only in A: {len(only_a)},  only in B: {len(only_b)}")
    if sa == sb:
        verdict = "identical (A = B)"
    elif sa < sb:
        verdict = "A is a subset of B"
    elif sb < sa:
        verdict = "B is a subset of A"
    else:
        verdict = "neither is a subset of the other"
    print(f"  -> {verdict}")
    for label, keys, src in (("only in A", only_a, a), ("only in B", only_b, b)):
        if keys:
            print(f"  {label} ({len(keys)}):   [ID  rest of the line]")
            for k in keys:
                print(f"    {k}  {src[k]}")
    print()
    return verdict


def main():
    args = sys.argv[1:]
    if len(args) != 2:
        print(__doc__)
        sys.exit(1)
    a, b = args
    if os.path.isfile(a) and os.path.isfile(b):                 # two files of the same kind
        what = "bad runs" if "badruns" in a else "bad towers" if "badtowers" in a else "entries"
        compare(a, b, what)
    else:                                                        # two outTags
        for suffix, what in (("_final_badtowers.txt", "final bad towers"), ("_final_badruns.txt", "final bad runs")):
            fa, fb = a + suffix, b + suffix
            missing = [p for p in (fa, fb) if not os.path.isfile(p)]
            if missing:
                print(f"=== {what} ===\n  missing: {', '.join(missing)}\n")
                continue
            compare(fa, fb, what)


if __name__ == "__main__":
    main()
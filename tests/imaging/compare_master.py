#!/usr/bin/env python3
"""Compare two ILLUMINA result sets pointing by pointing.

A result set is a file of ``key=value`` blocks separated by blank lines
(``<basenm>_results.txt`` written by the kernel, ``sky_results.txt`` of the
imaging scenario, or ``results.txt`` written by ``run_pointings.py``).
Every block must hold ``elevation_deg`` and ``azimuth_deg``; the pointings
are matched on these two values (rounded to 0.01 deg). Blocks that occur
several times (the grid holds azimuth 180 twice) are compared once.

For every compared quantity the script prints the number of pointings
where both values are zero, the maximum and median relative difference
``|b - a| / max(|a|, |b|)`` over the pointings where at least one value is
non-zero, and the pointings with the largest differences. ``los_steps``
(written by ``run_pointings.py``) is shown when both sets have it.

    python3 tests/imaging/compare_master.py A/results.txt B/results.txt \\
        --label-a master --label-b ai-update --top 12 [--markdown out.md]
"""
import argparse
import sys

DEFAULT_KEYS = ["diffuse_radiance", "direct_radiance_sources"]
ALL_KEYS = [
    "direct_irradiance_sources",
    "direct_irradiance_reflection",
    "direct_radiance_sources",
    "direct_radiance_reflection",
    "cloud_radiance",
    "diffuse_radiance",
]


def read_results(path):
    """Return {(elev, azim): {key: value}}."""
    out = {}
    block = {}
    with open(path) as f:
        lines = list(f) + [""]
    for ln in lines:
        ln = ln.strip()
        if not ln:
            if block:
                key = (round(block["elevation_deg"], 2), round(block["azimuth_deg"], 2))
                out.setdefault(key, block)
                block = {}
            continue
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        try:
            block[k.strip()] = float(v)
        except ValueError:
            block[k.strip()] = v.strip()
    return out


def rel_diff(a, b):
    m = max(abs(a), abs(b))
    return 0.0 if m == 0 else abs(b - a) / m


def compare(ra, rb, keys):
    common = sorted(set(ra) & set(rb))
    rows = []
    for p in common:
        row = {"pointing": p}
        for k in keys:
            a, b = ra[p].get(k), rb[p].get(k)
            row[k] = (a, b, None if a is None or b is None else rel_diff(a, b))
        if "los_steps" in ra[p] and "los_steps" in rb[p]:
            row["los_steps"] = (int(ra[p]["los_steps"]), int(rb[p]["los_steps"]))
        rows.append(row)
    return common, rows


def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def summarise(rows, keys):
    summary = {}
    for k in keys:
        diffs = [(r[k][2], r) for r in rows if r[k][2] is not None and (r[k][0] != 0 or r[k][1] != 0)]
        both_zero = sum(1 for r in rows if r[k][2] is not None and r[k][0] == 0 and r[k][1] == 0)
        d = [x for x, _ in diffs]
        summary[k] = {
            "n": len(rows),
            "both_zero": both_zero,
            "n_diff": sum(1 for x in d if x > 0),
            "max": max(d) if d else 0.0,
            "median": median(d),
            "top": sorted(diffs, key=lambda t: -t[0]),
        }
    return summary


def fmt_pointing(p):
    return "el %5.1f az %6.1f" % p


def print_report(la, lb, na, nb, common, rows, summary, keys, top, md=None):
    out = []
    out.append("%d pointings in %s, %d in %s, %d compared" % (na, la, nb, lb, len(common)))
    out.append("")
    out.append("| quantity | compared | both zero | differing | max rel diff | median rel diff |")
    out.append("|----------|---------:|----------:|----------:|-------------:|----------------:|")
    for k in keys:
        s = summary[k]
        out.append("| %s | %d | %d | %d | %.2e | %.2e |" % (k, s["n"], s["both_zero"], s["n_diff"], s["max"], s["median"]))
    steps = [r for r in rows if "los_steps" in r and r["los_steps"][0] != r["los_steps"][1]]
    if any("los_steps" in r for r in rows):
        out.append("")
        out.append("line of sight steps differ on %d of %d pointings" % (len(steps), len(rows)))
    for k in keys:
        out.append("")
        out.append("Largest relative differences of %s (%s -> %s):" % (k, la, lb))
        out.append("")
        head = "| pointing | %s | %s | rel diff |" % (la, lb)
        sep = "|----------|---:|---:|---:|"
        if any("los_steps" in r for r in rows):
            head += " steps %s | steps %s |" % (la, lb)
            sep += "---:|---:|"
        out.append(head)
        out.append(sep)
        for d, r in summary[k]["top"][:top]:
            a, b, _ = r[k]
            line = "| %s | %.3e | %.3e | %.2e |" % (fmt_pointing(r["pointing"]), a, b, d)
            if "los_steps" in r:
                line += " %d | %d |" % r["los_steps"]
            elif any("los_steps" in x for x in rows):
                line += " | |"
            out.append(line)
    text = "\n".join(out)
    print(text)
    if md:
        with open(md, "w") as f:
            f.write(text + "\n")
        print("\nwritten:", md)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", help="first result set (reference, e.g. master)")
    ap.add_argument("b", help="second result set (e.g. ai-update)")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--keys", default=",".join(DEFAULT_KEYS),
                    help="comma separated quantities (default %s; 'all' for the six)" % ",".join(DEFAULT_KEYS))
    ap.add_argument("--top", type=int, default=10, help="pointings listed per quantity")
    ap.add_argument("--markdown", default=None, help="also write the report to this file")
    ap.add_argument("--rtol", type=float, default=None,
                    help="exit with status 1 when a relative difference exceeds this value")
    args = ap.parse_args(argv)

    keys = ALL_KEYS if args.keys == "all" else args.keys.split(",")
    ra, rb = read_results(args.a), read_results(args.b)
    common, rows = compare(ra, rb, keys)
    if not common:
        sys.exit("no common pointing between the two result sets")
    summary = summarise(rows, keys)
    print_report(args.label_a, args.label_b, len(ra), len(rb), common, rows, summary, keys, args.top, args.markdown)
    if args.rtol is not None and any(summary[k]["max"] > args.rtol for k in keys):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

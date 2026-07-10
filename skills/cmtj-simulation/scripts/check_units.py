#!/usr/bin/env python3
"""Heuristic static checker for common CMTJ unit/parameter mistakes.

Not a physics validator -- just greps a script for the mistakes that show up
most often in AI-generated cmtj code (see AGENTS.md unit table). Regex-based
on purpose: works on any .py without importing cmtj or running the sim.

Usage: python check_units.py <file.py> [file2.py ...]
Exit code: 1 if any warning was printed, else 0.
"""

import re
import sys

MS_RE = re.compile(r"\bMs\s*=\s*([0-9.eE+-]+)")
THICKNESS_RE = re.compile(r"\bthickness\s*=\s*([0-9.eE+-]+)")
DT_RE = re.compile(r"\bdt\s*=\s*([0-9.eE+-]+)")
J_RE = re.compile(r"\b(?:J1|J2|J)\s*=\s*([0-9.eE+-]+)")
LAYER_CTOR_RE = re.compile(r"\bLayer\s*\(")
FIXED_LAYER_NAME_RE = re.compile(r"""["'](reference|fixed|pinned)["']""", re.IGNORECASE)
SB_FILE_HINT = re.compile(r"general_sb")


def check(path: str) -> int:
    with open(path) as f:
        text = f.read()

    is_sb = bool(SB_FILE_HINT.search(text))
    warnings = []

    for m in MS_RE.finditer(text):
        val = float(m.group(1))
        if is_sb:
            if val < 1e4:
                warnings.append(
                    f"Ms={val}: looks like Tesla, but this looks like a "
                    "Smit-Beljers file (Ms should be A/m, ~1e5-1e6)"
                )
        else:
            if val > 10:
                warnings.append(f"Ms={val}: looks like A/m, but core Layer API wants Tesla (~0.5-1.6)")

    for m in THICKNESS_RE.finditer(text):
        val = float(m.group(1))
        if not (1e-11 <= val <= 1e-7):
            warnings.append(f"thickness={val}: outside typical 0.8-2.0 nm range (in meters, ~1e-9)")

    if FIXED_LAYER_NAME_RE.search(text) and "setReferenceLayer" not in text:
        warnings.append(
            "layer named reference/fixed/pinned but no setReferenceLayer() call -- "
            "fixed/reference layers should be modelled via layer.setReferenceLayer(...) on the "
            "free layer, not as a separate Layer object, see AGENTS.md"
        )

    dt_vals = [float(m.group(1)) for m in DT_RE.finditer(text)]
    j_vals = [float(m.group(1)) for m in J_RE.finditer(text)]
    if dt_vals and j_vals:
        max_dt = max(dt_vals)
        max_j = max(abs(v) for v in j_vals)
        if (
            max_j >= 1e-4
            and max_dt > 1e-13
            and "DormandPrince" not in text
            and "adaptive" not in text.lower()
        ):
            warnings.append(
                f"dt={max_dt} with IEC J={max_j}: large IEC usually needs "
                "dt <= 1e-13 for stable fixed-step integration"
            )

    for w in warnings:
        print(f"{path}: {w}")
    return len(warnings)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 0
    total = sum(check(p) for p in sys.argv[1:])
    return 1 if total else 0


if __name__ == "__main__":
    raise SystemExit(main())

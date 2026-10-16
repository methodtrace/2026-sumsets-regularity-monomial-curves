"""Run every check and write its output to results/<name>.txt.

    python3 run_all.py            # run all checks, write results/, compare with MANIFEST.sha256
    python3 run_all.py --write    # also rewrite MANIFEST.sha256 from the new results

Each check is deterministic: the files in results/ depend only on the code, not on the machine,
so their SHA-256 digests must equal those in MANIFEST.sha256. Timing goes to the terminal only.
"""
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHECKS = [
    "selftest",
    "theorem_a",
    "case3_identities",
    "theorem_b",
    "family_p3",
    "sigma_literal",
    "defect",
    "band_family",
    "gluing",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if sys.flags.optimize:
        sys.exit("run without -O: the checks rely on assert statements")
    env = dict(os.environ)
    env.pop("PYTHONOPTIMIZE", None)
    env["PYTHONIOENCODING"] = "utf-8"
    results = ROOT / "results"
    results.mkdir(exist_ok=True)
    lines = []
    for name in CHECKS:
        start = time.time()
        out = subprocess.run([sys.executable, str(ROOT / "checks" / f"{name}.py")],
                             capture_output=True, cwd=ROOT, env=env)
        if out.returncode != 0:
            sys.stderr.write(out.stderr.decode("utf-8", "replace"))
            sys.exit(f"FAILED: {name}")
        target = results / f"{name}.txt"
        target.write_bytes(out.stdout.replace(b"\r\n", b"\n"))
        lines.append(f"{digest(target)}  results/{name}.txt")
        print(f"{name:18s} ok  {time.time() - start:7.1f} s", file=sys.stderr)
    manifest = ROOT / "MANIFEST.sha256"
    new = "\n".join(lines) + "\n"
    if "--write" in sys.argv:
        manifest.write_bytes(new.encode("ascii"))
        print("MANIFEST.sha256 written", file=sys.stderr)
    elif manifest.exists():
        if manifest.read_bytes().replace(b"\r\n", b"\n") == new.encode("ascii"):
            print("all results match MANIFEST.sha256 byte for byte", file=sys.stderr)
        else:
            old = set(manifest.read_bytes().decode("ascii").splitlines())
            for line in lines:
                if line not in old:
                    print(f"MISMATCH: {line.split()[1]}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()

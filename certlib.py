"""Small helpers shared by the verifiers: argument handling, input hashes and
deterministic JSON output.  Standard library only."""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def input_hashes(script_name, extra=("mprk_exact.py", "mprk_maps.py", "certlib.py")):
    names = [script_name] + [e for e in extra if e != script_name]
    return {name: sha256_file(HERE / name) for name in names}


def parameters_hash(params):
    blob = json.dumps(params, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def out_path():
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--out":
        p = Path(args[1])
        if p.exists():
            raise SystemExit("refusing to overwrite existing file %s" % p)
        return p
    raise SystemExit("usage: python -S -O %s --out NEW_CERTIFICATE.json" % Path(sys.argv[0]).name)


def write_certificate(path, result):
    text = json.dumps(result, indent=2) + "\n"
    path.write_text(text)
    return text

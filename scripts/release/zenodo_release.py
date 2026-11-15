"""Archive a tagged eco_opt release on Zenodo as a new version of the record.

Usage (from the repository root, after bumping the version, committing,
tagging ``v<version>`` and pushing both the branch and the tag):

    python scripts/release/zenodo_release.py            # create, upload, publish
    python scripts/release/zenodo_release.py --draft    # stop before publishing

The Zenodo token is read from ``ZENODO_API_KEY`` or from
``~/.config/zenodo/api.env`` (line ``ZENODO_API_KEY=...``). It is never
written to the repository. The record state is kept in
``scripts/release/zenodo.json``; commit that file after a release.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

API = "https://zenodo.org/api"
ROOT = Path(__file__).resolve().parents[2]
STATE = ROOT / "scripts" / "release" / "zenodo.json"


def token() -> str:
    tok = os.environ.get("ZENODO_API_KEY")
    if tok:
        return tok.strip()
    env = Path.home() / ".config" / "zenodo" / "api.env"
    for line in env.read_text().splitlines():
        if line.startswith("ZENODO_API_KEY="):
            return line.split("=", 1)[1].strip()
    sys.exit("No Zenodo token found.")


def call(method: str, url: str, data=None, raw: bytes | None = None, ctype="application/json"):
    body = raw if raw is not None else (json.dumps(data).encode() if data is not None else None)
    req = urllib.request.Request(url, data=body, method=method,
                                 headers={"Authorization": "Bearer " + token(), "Content-Type": ctype})
    with urllib.request.urlopen(req) as r:
        text = r.read()
    return json.loads(text) if text else None


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()


def project_version() -> str:
    m = re.search(r'^version\s*=\s*"([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M)
    return m.group(1)


def checks(version: str) -> str:
    if git("status", "--porcelain"):
        sys.exit("Working tree is not clean.")
    cff = re.search(r'^version:\s*"?([^"\n]+)"?', (ROOT / "CITATION.cff").read_text(), re.M).group(1)
    if cff != version:
        sys.exit(f"CITATION.cff version {cff} differs from pyproject version {version}.")
    init = re.search(r'^__version__\s*=\s*"([^"]+)"', (ROOT / "src" / "eco_opt" / "__init__.py").read_text(), re.M).group(1)
    if init != version:
        sys.exit(f"eco_opt.__version__ {init} differs from pyproject version {version}.")
    tag = f"v{version}"
    if tag not in git("tag", "--list", tag).split():
        sys.exit(f"Tag {tag} does not exist.")
    git("fetch", "--tags", "origin")
    if git("rev-parse", f"{tag}^{{commit}}") != git("rev-parse", "HEAD"):
        sys.exit(f"HEAD is not at tag {tag}.")
    if git("rev-parse", "HEAD") != git("rev-parse", "origin/main"):
        sys.exit("HEAD is not pushed to origin/main.")
    if not git("ls-remote", "--tags", "origin", tag):
        sys.exit(f"Tag {tag} is not pushed.")
    return tag


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draft", action="store_true", help="upload but do not publish")
    ap.add_argument("--deposition", help="use this existing draft deposition id (first release)")
    args = ap.parse_args()

    state = json.loads(STATE.read_text())
    version = project_version()
    tag = checks(version)
    if version == state["latest_version"] and not args.deposition:
        sys.exit(f"Version {version} is already archived ({state['latest_doi']}). Bump the version first.")

    if args.deposition:
        dep = call("GET", f"{API}/deposit/depositions/{args.deposition}")
    else:
        nv = call("POST", f"{API}/deposit/depositions/{state['latest_deposition_id']}/actions/newversion")
        dep = call("GET", nv["links"]["latest_draft"])
        for f in dep.get("files", []):
            call("DELETE", f"{API}/deposit/depositions/{dep['id']}/files/{f['id']}")

    meta = dep["metadata"]
    meta["version"] = version
    meta["publication_date"] = dt.date.today().isoformat()
    meta.pop("doi", None)
    # The legacy API reports BSD-2-Clause as "bsd-2-clause-netbsd"; always send the SPDX id.
    meta["license"] = "bsd-2-clause"
    dep = call("PUT", f"{API}/deposit/depositions/{dep['id']}", {"metadata": meta})

    name = f"eco_opt-{version}.zip"
    with tempfile.TemporaryDirectory() as tmp:
        zpath = Path(tmp) / name
        subprocess.run(["git", "archive", "--format=zip", f"--prefix=eco_opt-{version}/",
                        "-o", str(zpath), tag], cwd=ROOT, check=True)
        call("PUT", f"{dep['links']['bucket']}/{name}", raw=zpath.read_bytes(),
             ctype="application/octet-stream")

    doi = dep["metadata"]["prereserve_doi"]["doi"]
    if args.draft:
        print(f"Draft {dep['id']} ready (DOI {doi} reserved): {dep['links']['html']}")
        return
    pub = call("POST", f"{API}/deposit/depositions/{dep['id']}/actions/publish")
    state.update(latest_deposition_id=str(pub["id"]), latest_version=version, latest_doi=pub["doi"])
    STATE.write_text(json.dumps(state, indent=2) + "\n")
    print(f"Published {tag}: https://doi.org/{pub['doi']} (all versions: https://doi.org/{state['concept_doi']})")
    print("Commit scripts/release/zenodo.json, and update the DOI in the manuscript if it cites this version.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Resumable, checksum-verifying fetcher for the VCC 2026 data plan.

Reads plan rows from manifests/fetch_plan.csv, downloads with .partial +
byte-range resume, verifies size and (when the source published one) the
checksum, and records every attempt in manifests/acquisition_manifest.csv.

Nothing here invents a checksum: `checksum_source` is empty when the origin
published none, and in that case `verified` stays 'size_only'.
"""
import argparse, csv, hashlib, json, os, pathlib, subprocess, sys, time
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parent.parent
PLAN = ROOT / "manifests" / "fetch_plan.csv"
ACQ = ROOT / "manifests" / "acquisition_manifest.csv"
LOGDIR = ROOT / "logs"
ACQ_FIELDS = ["file_id","dataset_id","filename","role","source_url","release","expected_bytes",
              "actual_bytes","checksum_source","checksum_local","verified","local_path",
              "acquisition_status","attempts","last_error","updated_utc","resume_hint"]

def now(): return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def read_csv(p):
    if not p.exists(): return []
    with open(p, encoding="utf-8-sig", newline="") as f: return list(csv.DictReader(f))

def write_csv(p, rows, fields):
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    os.replace(tmp, p)

def digest(path, algo, chunk=16 << 20):
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        while (b := f.read(chunk)): h.update(b)
    return h.hexdigest().lower()

def looks_like_html(path):
    """Guard against a login/error page saved under a data filename."""
    with open(path, "rb") as f: head = f.read(512).lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html") or b"<title>" in head[:200]

def curl(url, dest, expected, log, connections=1, timeout=0):
    """Download url -> dest.partial, resuming when the file already exists."""
    part = dest.with_name(dest.name + ".partial")
    start = part.stat().st_size if part.exists() else 0
    if expected and start > expected:
        return False, f"partial larger than expected ({start} > {expected})"
    if expected and start == expected:
        os.replace(part, dest); return True, ""
    cmd = ["curl","-sS","-L","--fail","--retry","5","--retry-delay","5","--retry-all-errors",
           "--connect-timeout","30","--speed-time","120","--speed-limit","1024",
           "-C","-","-o",str(part),url]
    if timeout: cmd[1:1] = ["-m", str(timeout)]
    with open(log, "ab") as lf:
        lf.write(f"\n[{now()}] GET {url} -> {part} (resume from {start})\n".encode())
        r = subprocess.run(cmd, stdout=lf, stderr=lf)
    if r.returncode != 0:
        return False, f"curl exit {r.returncode}"
    got = part.stat().st_size if part.exists() else 0
    if expected and got != expected:
        return False, f"size mismatch: got {got}, expected {expected}"
    if got < 4096 or looks_like_html(part):
        if looks_like_html(part):
            return False, "response body is HTML (login or error page), not data"
    os.replace(part, dest)
    return True, ""

def one(row, acq_by_id, force=False):
    fid = row["file_id"]; dest = ROOT / row["local_path"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    exp = int(row["expected_bytes"]) if row.get("expected_bytes") else 0
    csum = (row.get("checksum_source") or "").strip()
    rec = acq_by_id.get(fid, {})
    attempts = int(rec.get("attempts") or 0)
    log = LOGDIR / f"fetch_{row['dataset_id']}.log"; LOGDIR.mkdir(exist_ok=True)

    def finish(status, verified, err="", local_hash=""):
        out = dict(row); out.update(
            actual_bytes=dest.stat().st_size if dest.exists() else 0,
            checksum_local=local_hash, verified=verified, acquisition_status=status,
            attempts=attempts, last_error=err, updated_utc=now(),
            resume_hint=f"python scripts/fetch.py --ids {fid}")
        acq_by_id[fid] = out
        return out

    # Already complete and verified? Skip without re-hashing a large file.
    prev = rec.get("acquisition_status")
    if dest.exists() and not force:
        size_ok = (not exp) or dest.stat().st_size == exp
        if size_ok and prev == "COMPLETE" and rec.get("verified") in ("checksum_match","size_only"):
            print(f"  skip {fid} {dest.name} (already {rec['verified']})"); return finish("COMPLETE", rec["verified"], "", rec.get("checksum_local",""))
        if not size_ok:
            print(f"  {fid}: existing file wrong size -> re-fetch")
            dest.rename(dest.with_name(dest.name + ".partial"))

    attempts += 1
    ok, err = curl(row["source_url"], dest, exp, log, timeout=int(row.get("timeout") or 0))
    if not ok:
        print(f"  FAIL {fid}: {err}")
        return finish("RETRYABLE_FAILURE", "unverified", err)
    lh = ""; verified = "size_only"
    if csum:
        algo, expect = csum.split(":", 1)
        lh = digest(dest, algo)
        if lh != expect.strip().lower():
            print(f"  CHECKSUM MISMATCH {fid}: {lh} != {expect}")
            return finish("RETRYABLE_FAILURE", "checksum_mismatch", "checksum mismatch", lh)
        verified = "checksum_match"
    else:
        lh = "sha256:" + digest(dest, "sha256")
    print(f"  ok   {fid} {dest.name} {dest.stat().st_size} {verified}")
    return finish("COMPLETE", verified, "", lh)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", help="comma-separated file_id list")
    ap.add_argument("--dataset", help="fetch every planned file of this dataset_id")
    ap.add_argument("--group", help="fetch every planned file of this fetch_group")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    plan = read_csv(PLAN)
    sel = plan
    if a.ids: want = set(a.ids.split(",")); sel = [r for r in plan if r["file_id"] in want]
    elif a.dataset: sel = [r for r in plan if r["dataset_id"] == a.dataset]
    elif a.group: sel = [r for r in plan if r.get("fetch_group") == a.group]
    if not sel: sys.exit("no planned files matched")
    total = sum(int(r["expected_bytes"] or 0) for r in sel)
    print(f"selected {len(sel)} files, {total/1e9:.3f} GB expected")
    if a.dry_run:
        for r in sel: print(f"  {r['file_id']:6s} {int(r['expected_bytes'] or 0)/1e9:9.3f} GB  {r['local_path']}")
        return
    acq = {r["file_id"]: r for r in read_csv(ACQ)}
    for r in sel:
        one(r, acq, force=a.force)
        write_csv(ACQ, list(acq.values()), ACQ_FIELDS)
    done = sum(1 for r in acq.values() if r["acquisition_status"] == "COMPLETE")
    print(f"manifest: {done}/{len(acq)} COMPLETE -> {ACQ}")

if __name__ == "__main__": main()

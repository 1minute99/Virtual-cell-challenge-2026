"""Restore the 157 frozen Tahoe-100M parquet shards used by the context_signal / drug_context_pretrain pipeline.

Source of the exact list: vcc2026_data/reports/context_signal_2026-09-15/tahoe/candidate_shard_plan.json
(157 entries: path, bytes, sha256), revision 2dc57900b7981cfcf5e211527169a0b006546a95. Each entry's sha256
and size are cross-checked against the HF LFS listing saved at
reports/context_signal_2026-09-15/source_documents/tahoe_data_directory.json before anything is downloaded.

Downloads https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/<rev>/<path> into
reports/context_signal_2026-09-15/tahoe/raw/<path>. No credentials. Resumable (HTTP Range into <file>.part),
skips files whose sha256 already matches, verifies every file's sha256 and size, and exits non-zero on any
mismatch (a mismatching .part is deleted so the next run starts that file fresh).
Writes a receipt JSON to reports/context_signal_2026-09-15/tahoe/raw/restore_receipt_<date>.json.
"""
from __future__ import annotations
import hashlib, json, os, sys, time, urllib.request, urllib.error
from datetime import datetime, timezone
from pathlib import Path

REV = '2dc57900b7981cfcf5e211527169a0b006546a95'
BASE = Path('/home/wonmin-kim/AIVC/vcc2026_data/reports/context_signal_2026-09-15')
PLAN = BASE / 'tahoe/candidate_shard_plan.json'
LISTING = BASE / 'source_documents/tahoe_data_directory.json'
DEST = BASE / 'tahoe/raw'
URL = 'https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/{rev}/{path}'
CHUNK = 8 << 20


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(CHUNK), b''):
            h.update(b)
    return h.hexdigest()


def log(**kw):
    kw['t'] = datetime.now(timezone.utc).strftime('%H:%M:%S')
    print(json.dumps(kw), flush=True)


def load_manifest():
    plan = json.loads(PLAN.read_text())
    assert plan['revision'] == REV, 'plan revision differs'
    shards = plan['shards']
    assert len(shards) == 157, len(shards)
    listing = {x['path']: x for x in json.loads(LISTING.read_text())}
    for s in shards:
        L = listing[s['path']]
        if L['lfs']['oid'] != s['sha256'] or int(L['size']) != int(s['bytes']):
            raise SystemExit(f"manifest disagreement for {s['path']}: plan vs HF listing")
    return shards


def download(s, retries=6):
    out = DEST / s['path']
    out.parent.mkdir(parents=True, exist_ok=True)
    part = out.with_name(out.name + '.part')
    url = URL.format(rev=REV, path=s['path'])
    size = int(s['bytes'])
    for attempt in range(1, retries + 1):
        have = part.stat().st_size if part.exists() else 0
        if have > size:
            part.unlink(); have = 0
        if have < size:
            req = urllib.request.Request(url, headers={'User-Agent': 'aivc-restore/1.0'})
            if have:
                req.add_header('Range', f'bytes={have}-')
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    if have and r.status != 206:      # server ignored Range: start over
                        have = 0
                        mode = 'wb'
                    else:
                        mode = 'ab' if have else 'wb'
                    with open(part, mode) as f:
                        while True:
                            b = r.read(CHUNK)
                            if not b:
                                break
                            f.write(b)
            except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
                log(event='retry', path=s['path'], attempt=attempt, err=repr(e)[:200])
                time.sleep(min(60, 5 * attempt))
                continue
        if part.stat().st_size != size:
            log(event='retry', path=s['path'], attempt=attempt, err='short file', got=part.stat().st_size)
            time.sleep(5)
            continue
        got = sha256_file(part)
        if got != s['sha256']:
            part.unlink()
            raise SystemExit(f"SHA256 MISMATCH {s['path']}: expected {s['sha256']} got {got}")
        os.replace(part, out)
        return 'downloaded'
    raise SystemExit(f"FAILED after {retries} attempts: {s['path']}")


def main():
    shards = load_manifest()
    t0 = time.time()
    rows = []
    for i, s in enumerate(shards, 1):
        out = DEST / s['path']
        status = None
        if out.exists() and out.stat().st_size == int(s['bytes']):
            if sha256_file(out) == s['sha256']:
                status = 'already_ok'
            else:
                log(event='bad_existing', path=s['path']); out.unlink()
        elif out.exists():
            log(event='bad_existing_size', path=s['path']); out.unlink()
        if status is None:
            status = download(s)
        # final verification of the file in place
        if out.stat().st_size != int(s['bytes']) or sha256_file(out) != s['sha256']:
            raise SystemExit(f"POST-VERIFY MISMATCH {s['path']}")
        rows.append(dict(path=s['path'], bytes=int(s['bytes']), sha256=s['sha256'], status=status, verified=True))
        log(event='shard', i=i, n=len(shards), path=s['path'], status=status,
            elapsed_s=round(time.time() - t0))
    receipt = dict(revision=REV, manifest=str(PLAN), listing=str(LISTING), dest=str(DEST),
                   n_expected=len(shards), n_verified=sum(r['verified'] for r in rows),
                   bytes_total=sum(r['bytes'] for r in rows), elapsed_s=round(time.time() - t0),
                   finished_utc=datetime.now(timezone.utc).isoformat(), shards=rows)
    rp = DEST / f"restore_receipt_{datetime.now().strftime('%Y-%m-%d')}.json"
    rp.write_text(json.dumps(receipt, indent=1))
    log(event='done', verified=receipt['n_verified'], expected=len(shards), receipt=str(rp))
    if receipt['n_verified'] != len(shards):
        sys.exit(1)


if __name__ == '__main__':
    main()

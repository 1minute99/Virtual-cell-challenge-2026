"""Optional exact-count accelerator; original scorer and receipts stay intact."""
from pathlib import Path
import json


def install_acceleration(scorer, accelerated, extra_hashes):
    original_worker=scorer.ContextForcingScoreWorker
    original_bulk=scorer.E._reconstructed_bulk
    original_hashes=scorer.implementation_hashes
    marker='context_forcing_fast_bulk.py'
    if marker not in extra_hashes:
        raise ValueError('Accelerator identity missing')

    class AcceleratedWorker(original_worker):
        def score(self,request):
            output=Path(request['output_path'])
            fast=True
            if output.exists():
                saved=json.loads(output.read_text())
                # Previously completed original receipts are verified by their
                # original implementation identity and returned without writes.
                fast=marker in saved.get('implementation_hashes',{})
            try:
                scorer.E._reconstructed_bulk=accelerated if fast else original_bulk
                scorer.implementation_hashes=(lambda:{**original_hashes(),**extra_hashes}) if fast else original_hashes
                return super().score(request)
            finally:
                scorer.E._reconstructed_bulk=original_bulk
                scorer.implementation_hashes=original_hashes

    scorer.ContextForcingScoreWorker=AcceleratedWorker


def main():
    # The legacy evaluator deliberately removes the scripts directory during
    # its imports; restore it before loading the other standalone entry point.
    import sys
    import context_forcing_fast_bulk as F
    here=str(Path(__file__).resolve().parent)
    if here not in sys.path:sys.path.insert(0,here)
    import context_forcing_score as S
    paths=[Path(F.__file__),Path(__file__)]
    hashes={p.name:S.sha256(p) for p in paths}
    install_acceleration(S,F.reconstructed_bulk,hashes)
    return S.main()


if __name__=='__main__':
    raise SystemExit(main())

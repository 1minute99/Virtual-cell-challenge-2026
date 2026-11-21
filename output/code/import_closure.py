"""Static import closure of VCC scripts (local modules only) + subprocess-called scripts named in source text.
Usage: python import_closure.py ROOT_SCRIPT [...]  -> prints one path per line (absolute)."""
import ast, pathlib, re, sys
SD = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
EXTRA_DIRS = [SD]
def local(mod):
    p = SD / (mod.split(".")[0] + ".py")
    return p if p.exists() else None
seen = set(); stack = [pathlib.Path(a).resolve() for a in sys.argv[1:]]
while stack:
    f = stack.pop()
    if f in seen or not f.exists(): continue
    seen.add(f)
    txt = f.read_text(errors="ignore")
    for m in re.finditer(r"scripts/([A-Za-z0-9_]+\.(?:py|sh))", txt):     # subprocess / shell references
        q = SD / m.group(1)
        if q.exists(): stack.append(q)
    if f.suffix != ".py": continue
    try: tree = ast.parse(txt)
    except SyntaxError: continue
    for n in ast.walk(tree):
        names = []
        if isinstance(n, ast.Import): names = [a.name for a in n.names]
        elif isinstance(n, ast.ImportFrom) and n.module and n.level == 0: names = [n.module]
        for nm in names:
            q = local(nm)
            if q: stack.append(q)
for f in sorted(seen): print(f)

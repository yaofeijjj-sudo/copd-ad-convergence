from pathlib import Path
import json

params = {}
for line in (Path.cwd() / "config/params.yml").read_text().splitlines():
    if not line.strip() or line.lstrip().startswith("#"):
        continue
    key, value = line.split(":", 1)
    try:
        params[key.strip()] = json.loads(value.strip())
    except json.JSONDecodeError:
        params[key.strip()] = value.strip()

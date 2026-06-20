from pathlib import Path

def get_project_root():
    path = Path(__file__).resolve()
    for parent in [path] + list(path.parents):
        if (parent / "pyproject.toml").exists():
            return parent
    raise RuntimeError("Project root not found")

PATH_ROOT = get_project_root()

DATA = PATH_ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
OUTPUTS = PATH_ROOT / "outputs"
from pathlib import Path
from datetime import datetime
import re
import shutil
import subprocess

project = Path(__file__).resolve().parent.parent
stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
run_dir = project / "repro_checks" / stamp

inputs = [
    "data/raw/GSE316391_counts_PE.csv.gz",
    "data/metadata/GSE316391_clinical_metadata.xlsx",
    "data/metadata/msigdb_hallmark_human.rds",
]

for name in inputs:
    if not (project / name).is_file():
        raise SystemExit(f"Missing source input: {name}")

numbered = sorted(
    path for path in (project / "scripts").glob("*.R")
    if re.match(r"^\d{2}_", path.name)
)

numbers = [int(path.name[:2]) for path in numbered]
if numbers != list(range(22)):
    raise SystemExit(
        f"Expected exactly one script for stages 00–21; found: {numbers}"
    )

for directory in [
    "scripts", "data/raw", "data/metadata", "data/processed",
    "results/tables", "results/figures", "docs", "logs"
]:
    (run_dir / directory).mkdir(parents=True, exist_ok=True)

for path in (project / "scripts").glob("*.R"):
    shutil.copy2(path, run_dir / "scripts" / path.name)

for name in inputs:
    shutil.copy2(project / name, run_dir / name)

# Plans are documentation inputs, not previous analysis outputs.
for name in ["docs/analysis_plan.md", "docs/qc_decisions.md"]:
    if (project / name).is_file():
        shutil.copy2(project / name, run_dir / name)

(project / "repro_checks" / "latest_run.txt").write_text(
    str(run_dir) + "\n", encoding="utf-8"
)

jobs = [
    ("check_environment.R", []),
    ("00_prepare_inputs.R", ["data/processed"]),
]
jobs += [(path.name, []) for path in numbered if path.name[:2] != "00"]

print(f"New run directory: {run_dir}", flush=True)

for filename, arguments in jobs:
    log = run_dir / "logs" / (Path(filename).stem + ".log")
    print(f"\nSTART: {filename}", flush=True)

    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            ["Rscript", f"scripts/{filename}", *arguments],
            cwd=run_dir,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )

    if result.returncode != 0:
        print(f"\nFAILED: {filename}", flush=True)
        print(f"Log: {log}", flush=True)
        lines = log.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
        print("\n".join(lines[-40:]), flush=True)
        raise SystemExit(result.returncode)

    print(f"COMPLETED: {filename}", flush=True)

(run_dir / "execution_completed.txt").write_text(
    "All analysis scripts completed successfully.\n",
    encoding="utf-8",
)

print("\nALL ANALYSIS SCRIPTS COMPLETED.", flush=True)
print("Numerical comparison with the original run is the next step.",
      flush=True)
print(f"Run directory: {run_dir}", flush=True)

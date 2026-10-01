from pathlib import Path
import shutil
import subprocess

project = Path(__file__).resolve().parent.parent
run_dir = Path(
    (project / "repro_checks/latest_run.txt").read_text().strip()
)

if not run_dir.is_dir():
    raise SystemExit("Yeniden çalıştırma dizini bulunamadı.")

scripts = [
    "16_leave_one_patient_out_models.R",
    "17_leave_one_patient_out_gsea.R",
    "18_hallmark_leave_one_out_heatmap.R",
]

source = project / "scripts" / scripts[0]
if "# Reset inherited normalization" not in source.read_text():
    raise SystemExit("Normalizasyon düzeltmesi scriptte bulunamadı.")

marker = run_dir / "execution_completed.txt"
if marker.exists():
    marker.unlink()

for name in scripts:
    shutil.copy2(
        project / "scripts" / name,
        run_dir / "scripts" / name
    )

    log = run_dir / "logs" / (
        Path(name).stem + "_normalization_fixed.log"
    )
    print(f"START: {name}", flush=True)

    with log.open("w") as handle:
        result = subprocess.run(
            ["Rscript", f"scripts/{name}"],
            cwd=run_dir,
            stdout=handle,
            stderr=subprocess.STDOUT
        )

    if result.returncode:
        print(f"FAILED: {name}\nLog: {log}", flush=True)
        print("\n".join(
            log.read_text(errors="replace").splitlines()[-40:]
        ))
        raise SystemExit(result.returncode)

    print(f"COMPLETED: {name}", flush=True)

marker.write_text(
    "Original full workflow completed; corrected stages 16-18 "
    "subsequently rerun successfully.\n"
)

print("CORRECTED SENSITIVITY RERUN COMPLETED.", flush=True)

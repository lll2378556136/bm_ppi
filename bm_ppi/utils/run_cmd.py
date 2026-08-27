from pathlib import Path
import subprocess


def run_cmd(cmd: str, outputs: list[Path], force: bool = False):
    for o in outputs:
        if o.exists() and not force:
            print(f"⏭️ {o} exists, skip")
            return
        o.parent.mkdir(parents=True, exist_ok=True)

    subprocess.run(cmd, shell=True, check=True)

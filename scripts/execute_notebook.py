"""Install a kernel inside .venv and execute the notebook in a fresh kernel."""
import os
from pathlib import Path
import subprocess
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root = Path(__file__).resolve().parents[1]
if Path(sys.prefix).resolve() != (root / ".venv").resolve():
    raise SystemExit("Execute with this project's .venv interpreter")
subprocess.run([sys.executable, "-m", "ipykernel", "install", "--prefix", str(root/".venv"),
                "--name", "q3hack-mvp", "--display-name", "Q3hack (.venv)"], check=True)
kernel_dir = root / ".venv" / "share" / "jupyter" / "kernels"
ksm = KernelSpecManager(kernel_dirs=[str(kernel_dir)])
assert Path(ksm.get_kernel_spec("q3hack-mvp").argv[0]).resolve() == Path(sys.executable).resolve()
runtime = root / "results" / "jupyter-runtime"
runtime.mkdir(parents=True, exist_ok=True)
os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime)
os.environ["JUPYTER_PATH"] = str(root / ".venv" / "share" / "jupyter")
km = KernelManager(kernel_name="q3hack-mvp", kernel_spec_manager=ksm)
path = root / "notebooks" / "mvp.ipynb"
nb = nbformat.read(path, as_version=4)
client = NotebookClient(nb, km=km, timeout=600, resources={"metadata": {"path": str(root)}},
                        allow_errors=False, record_timing=True)
client.execute()
nbformat.write(nb, path)
print(f"Executed {path}; all code cells succeeded in {sys.executable}")

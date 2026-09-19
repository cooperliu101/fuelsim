"""Compare native residuals and analytic tangents with independent PyTorch autograd."""

import argparse
from pathlib import Path
import subprocess
import sys
import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from model import load

p = argparse.ArgumentParser()
p.add_argument("--infer", required=True)
p.add_argument("--model", required=True)
p.add_argument("--work", type=Path, required=True)
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
torch.set_num_threads(1)
model, nodes, volume, signature = load(a.model)
rng = np.random.default_rng(919)
n = len(nodes)
x = model.input_mean.numpy() + rng.normal(0, 0.7, (12, n + 1)) * model.input_scale.numpy()
path = a.work / "inputs.txt"
np.savetxt(path, x, header=f"{len(x)} {n}", comments="", fmt="%.17g")
result = subprocess.run([a.infer, a.model, str(path)], check=True, text=True, capture_output=True)
values = np.array([[float(v) for v in row.split()] for row in result.stdout.splitlines()])
r = values[:, :n]
j = values[:, n:].reshape(-1, n, n)
reference = model(torch.tensor(x)).detach().numpy()
jref = np.array(
    [torch.autograd.functional.jacobian(model, torch.tensor(v)).detach().numpy()[:, :n] for v in x]
)
re = float(abs(r - reference).max())
je = float(abs(j - jref).max())
# C++ parser must reject truncated coefficients and unknown activations.
for name, text in [
    ("truncated", Path(a.model).read_text()[:500]),
    ("activation", Path(a.model).read_text().replace("tanh", "relu")),
]:
    bad = a.work / (name + ".txt")
    bad.write_text(text)
    if subprocess.run([a.infer, str(bad), str(path)], capture_output=True).returncode == 0:
        raise AssertionError("Malformed model accepted")
print(f"Python_CPP_residual_max_error_W={re:.12g} Python_CPP_tangent_max_error_W_K={je:.12g}")
if re > 1e-10 or je > 1e-12:
    raise AssertionError("Python/C++ mismatch")

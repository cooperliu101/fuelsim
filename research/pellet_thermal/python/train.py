"""CPU float64 tanh MLP, trained only against full FEM boundary residuals."""

import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
import torch
from model import PelletMLP, dataset, save


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--model", type=Path, required=True)
    args = p.parse_args()
    torch.set_num_threads(1)
    torch.manual_seed(1709)
    nodes, volume, signature, x, y = dataset(args.data)
    # Split before normalization; all four smooth-field families occur in each set.
    rng = np.random.default_rng(214)
    order = rng.permutation(len(x))
    train = order[: int(0.8 * len(x))]
    test = order[int(0.8 * len(x)) :]
    xm, xs = x[train].mean(0), x[train].std(0)
    ym, ys = y[train].mean(0), y[train].std(0)
    if np.any(xs <= 0) or np.any(ys <= 0):
        raise ValueError("Training set must excite every input and output")
    model = PelletMLP([x.shape[1], 64, y.shape[1]], xm, xs, ym, ys)
    tx, ty = torch.tensor(x[train]), torch.tensor(y[train])
    with torch.no_grad():
        # This first qualification is constant-k. Start tanh near its linear region,
        # without using the exact Schur matrix or its derivatives as training labels.
        model.layers[0].weight.normal_(0, 0.0002)
        model.layers[0].bias.zero_()
        features = torch.tanh(model.layers[0]((tx - model.input_mean) / model.input_scale))
        design = torch.cat([features, torch.ones((len(train), 1), dtype=torch.float64)], 1)
        fit = torch.linalg.lstsq(
            design, (ty - model.output_mean) / model.output_scale, rcond=1e-8, driver="gelsd"
        ).solution
        model.layers[1].weight.copy_(fit[:-1].T)
        model.layers[1].bias.copy_(fit[-1])
    optimizer = torch.optim.LBFGS(
        model.parameters(),
        max_iter=100,
        tolerance_grad=1e-12,
        tolerance_change=1e-15,
        line_search_fn="strong_wolfe",
    )

    def closure():
        optimizer.zero_grad()
        prediction = model(tx)
        mse = ((prediction - ty) / model.output_scale).square().mean()
        energy = ((prediction.sum(1) + tx[:, -1] * volume) / (tx[:, -1] * volume)).square().mean()
        loss = mse + 0.1 * energy
        loss.backward()
        return loss

    optimizer.step(closure)
    with torch.no_grad():
        prediction = model(torch.tensor(x[test])).numpy()
    difference = prediction - y[test]
    report = {
        "dataset_sha256": hashlib.sha256(args.data.read_bytes()).hexdigest(),
        "samples": len(x),
        "training_samples": len(train),
        "held_out_samples": len(test),
        "seed": 1709,
        "architecture": [x.shape[1], 64, y.shape[1]],
        "activation": "tanh",
        "dtype": "float64",
        "device": "cpu",
        "residual_relative_l2": float(np.linalg.norm(difference) / np.linalg.norm(y[test])),
        "residual_maximum_error_W": float(abs(difference).max()),
        "maximum_energy_error_W": float(abs(prediction.sum(1) + x[test, -1] * volume).max()),
        "torch_version": torch.__version__,
        "mesh_signature": str(signature),
        "training_temperature_min_K": float(x[train, :-1].min()),
        "training_temperature_max_K": float(x[train, :-1].max()),
        "training_source_min_W_m3": float(x[train, -1].min()),
        "training_source_max_W_m3": float(x[train, -1].max()),
    }
    # Autograd checks a complete 26-by-26 temperature tangent at a held-out sample.
    point = torch.tensor(x[test[0]], requires_grad=True)
    jac = torch.autograd.functional.jacobian(model, point).detach().numpy()[:, :-1]
    report["maximum_tangent_column_energy_error_W_K"] = float(abs(jac.sum(0)).max())
    if report["residual_relative_l2"] > 1e-4 or report["maximum_energy_error_W"] > 1e-3:
        raise RuntimeError(report)
    args.model.parent.mkdir(parents=True, exist_ok=True)
    save(args.model, model, nodes, volume, signature)
    args.model.with_suffix(".training.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Strict text format shared with the native C++ inference implementation."""

from pathlib import Path
import numpy as np
import torch


class PelletMLP(torch.nn.Module):
    def __init__(self, sizes, input_mean, input_scale, output_mean, output_scale):
        super().__init__()
        self.layers = torch.nn.ModuleList(
            [torch.nn.Linear(a, b, dtype=torch.float64) for a, b in zip(sizes[:-1], sizes[1:])]
        )
        for name, value in [
            ("input_mean", input_mean),
            ("input_scale", input_scale),
            ("output_mean", output_mean),
            ("output_scale", output_scale),
        ]:
            self.register_buffer(name, torch.as_tensor(value, dtype=torch.float64))

    def forward(self, x):
        x = (x - self.input_mean) / self.input_scale
        for i, layer in enumerate(self.layers):
            x = layer(x)
            if i + 1 < len(self.layers):
                x = torch.tanh(x)
        return x * self.output_scale + self.output_mean


def dataset(path):
    with Path(path).open() as f:
        if f.readline().split() != ["FUELSIM_PELLET_DATA", "1"]:
            raise ValueError("Dataset version mismatch")
        n, count, volume, signature = f.readline().split()
        nodes = np.array([[float(v) for v in f.readline().split()] for _ in range(int(n))])
        values = np.loadtxt(f, ndmin=2)
    if values.shape != (int(count), 2 * int(n) + 1) or not np.isfinite(values).all():
        raise ValueError("Dataset dimensions or values invalid")
    return nodes, float(volume), int(signature), values[:, : int(n) + 1], values[:, int(n) + 1 :]


def save(path, model, nodes, volume, signature):
    def row(f, x):
        f.write(" ".join(format(float(v), ".17g") for v in np.asarray(x).reshape(-1)) + "\n")

    with Path(path).open("w") as f:
        f.write(f"FUELSIM_PELLET_MLP 1\n{len(nodes)} {len(model.layers)} {volume:.17g} {signature}\n")
        for n in nodes:
            f.write(str(int(n[0])) + " " + " ".join(format(v, ".17g") for v in n[1:]) + "\n")
        for name in ["input_mean", "input_scale", "output_mean", "output_scale"]:
            row(f, getattr(model, name).detach().numpy())
        for i, layer in enumerate(model.layers):
            f.write(
                f"{layer.in_features} {layer.out_features} "
                + ("linear" if i + 1 == len(model.layers) else "tanh")
                + "\n"
            )
            row(f, layer.weight.detach().numpy())
            row(f, layer.bias.detach().numpy())


def load(path):
    tokens = iter(Path(path).read_text().split())
    if (next(tokens), next(tokens)) != ("FUELSIM_PELLET_MLP", "1"):
        raise ValueError("Model version mismatch")
    n, layers = int(next(tokens)), int(next(tokens))
    volume = float(next(tokens))
    signature = int(next(tokens))
    nodes = np.array([[float(next(tokens)) for _ in range(4)] for _ in range(n)])

    def array(count):
        return np.array([float(next(tokens)) for _ in range(count)])

    normalization = [array(n + 1), array(n + 1), array(n), array(n)]
    sizes = [n + 1]
    coefficients = []
    for i in range(layers):
        a, b, activation = int(next(tokens)), int(next(tokens)), next(tokens)
        if a != sizes[-1] or activation != ("linear" if i + 1 == layers else "tanh"):
            raise ValueError("Layer mismatch")
        sizes.append(b)
        coefficients.append((array(a * b).reshape(b, a), array(b)))
    if list(tokens):
        raise ValueError("Trailing model data")
    model = PelletMLP(sizes, *normalization)
    with torch.no_grad():
        for layer, (w, b) in zip(model.layers, coefficients):
            layer.weight.copy_(torch.from_numpy(w))
            layer.bias.copy_(torch.from_numpy(b))
    return model, nodes, volume, signature

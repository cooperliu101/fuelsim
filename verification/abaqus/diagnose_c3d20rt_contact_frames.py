"""Verify the convected first tangent on all native finite-friction frames.

The initial interface is normal to global X, so the conventional first tangent
is global Z. On this mesh it coincides with the reference eta tangent. Average
both current unit tangents by current sample area, then form the contact frame.
"""
import numpy as np
import diagnose_c3d20rt_finite_friction as d

first = np.zeros((len(d.X), 13, 3))
second = np.zeros_like(first)
for ids in d.faces:
    face = d.X[:, ids]
    target = [d.lookup[int(node)] for node in ids]
    for i, row in enumerate(target):
        rule = d.samples(i)
        shapes = np.array([d.shape(2*a-1, 2*b-1) for a, b, weight in rule])
        tx = np.einsum('qn,tnc->tqc', shapes[:, 1], face)
        ty = np.einsum('qn,tnc->tqc', shapes[:, 2], face)
        normal = -np.cross(tx, ty)
        area = np.linalg.norm(normal, axis=2)
        normal /= area[:, :, None]
        weights = 4 * (1 if i < 4 else 5) / 24 * rule[:, 2]
        t1 = ty / np.linalg.norm(ty, axis=2)[:, :, None]
        t2 = np.cross(normal, t1)
        first[:, row] += np.einsum('q,tq,tqc->tc', weights, area, t1)
        second[:, row] += np.einsum('q,tq,tqc->tc', weights, area, t2)
first /= np.linalg.norm(first, axis=2)[:, :, None]
second -= np.sum(first * second, axis=2)[:, :, None] * first
second /= np.linalg.norm(second, axis=2)[:, :, None]
first_error = np.linalg.norm(first-d.basis[:, :, 0], axis=2).max()
normal_error = np.linalg.norm(np.cross(first, second)-d.normal, axis=2).max()
print('convected_first_maximum_vector_difference', first_error)
print('convected_normal_maximum_vector_difference', normal_error)
print('native_output_tangent_maximum_nonorthogonality',
      np.abs(np.sum(d.basis[:, :, 0]*d.basis[:, :, 1], axis=2)).max())
assert first_error < 1e-7 and normal_error < 1e-7

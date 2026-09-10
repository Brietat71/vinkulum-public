"""Configuration derivatives of fixed-world wrenches on admitted R/P trees.

Original Jacobian cross-product assembly; no library Hessian tensor is used.
Derivation, tensor axes and independent references: docs/PINOCCHIO_LOADS.md.
"""

import numpy as np


def ancestor_coordinates(links):
    """Coordinate columns along each root-to-body path, in traversal order."""
    paths = {None: ()}
    for i, link in enumerate(links):
        paths[link.child] = (*paths[link.parent], i)
    return paths


def world_load_derivative(point_jacobian, angular_jacobian, force, moment, path):
    """Return d(tau_external[i])/dq[j] with time and world wrench fixed.

    For ordered ancestors j <= i, d(J_point[:,i])/dq[j] = W[:,j] × V[:,i].
    The point Hessian is symmetric in i,j. In contrast, d(W[:,i])/dq[j]
    is W[:,j] × W[:,i] only for strict ancestors j < i; it is zero otherwise.
    Prismatic W columns vanish. Coordinates on other branches contribute zero.
    Point transport is already included in the point Jacobian V.
    """
    indices = np.asarray(path, dtype=int)
    linear = point_jacobian[:, indices].T
    angular = angular_jacobian[:, indices].T
    # Array axes below are [effort row i, differentiated coordinate j, xyz].
    lower_force = np.tril(np.cross(angular[None, :, :], linear[:, None, :]) @ force)
    lower_moment = np.tril(
        np.cross(angular[None, :, :], angular[:, None, :]) @ moment, -1
    )
    block = lower_force + np.tril(lower_force, -1).T + lower_moment
    n = point_jacobian.shape[1]
    derivative = np.zeros((n, n))
    derivative[np.ix_(indices, indices)] = block
    return derivative

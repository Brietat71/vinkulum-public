"""Relecture binary64 du post-traitement des archives de repères 0.14.1–0.18.0.

L'archive a été produite avec des produits 3×3 accumulés par FMA, des
produits matrice-vecteur non fusionnés et une norme à somme séquentielle.
NumPy/BLAS choisit ces opérations selon le CPU : rejouer `@` ne garantit
donc pas les mêmes bits. Ici chaque arrondi est explicite, même sans FMA
matérielle (math.fma). Ce protocole relit une observation historique ; il
ne change ni le solveur, ni les sources figées, ni les critères rationnels.
"""
from math import fma, sqrt


def _produit3(a, b):
    return fma(a[2], b[2], fma(a[1], b[1], a[0] * b[0]))


def _matmul(a, b):
    columns = list(zip(*b))
    return [[_produit3(row, column) for column in columns] for row in a]


def ecarts(a, b, monde, matiere):
    """Vecteurs rotation (9) et vitesse (3), avec les arrondis historiques."""
    world = list(zip(*monde))
    material = list(zip(*matiere))
    native = b[2][0]
    rotation = _matmul(_matmul(world, [native[i:i + 3] for i in (0, 3, 6)]), material)
    dr = [x - y for x, y in zip([x for row in rotation for x in row], a[2][0])]
    w = b[3][0]
    # Ne pas utiliser sum : Python peut compenser les sommes flottantes.
    dw = [(row[0] * w[0] + row[1] * w[1]) + row[2] * w[2] - y
          for row, y in zip(world, a[3][0])]
    return dr, dw


def norme(values):
    square = 0.0
    for value in values:
        square += value * value
    return sqrt(square)

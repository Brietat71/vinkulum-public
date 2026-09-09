"""Isole l'arrondi des positions sur le premier palier de Princeton.

Les rotations du résultat habituel restent fixées. La reconstruction ne
résout QUE les équations de translation de la poutre « milieu », pas les
moments ni le problème statique complet. Elle n'est pas réinjectée au noyau.
NumPy longdouble doit offrir plus de précision que float64 sur la machine.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
from vinkulum import _vinkulum

from confronte_mbdyn import princeton
from diagnostic_statut_statique import appelle

ROOT = Path(__file__).resolve().parents[1]
CN = (2.84191e6, 6.40131e5, 9.03881e5)


def exp_rotation(v):
    """Rodrigues ; même seuil de série que le noyau, calcul NumPy séparé."""
    dtype = v.dtype.type
    square = v @ v
    skew = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]], dtype=dtype)
    if square < dtype(1e-8):
        a = 1 - square / 6 + square * square / 120
        b = dtype(.5) - square / 24 + square * square / 720
    else:
        angle = np.sqrt(square)
        a = np.sin(angle) / angle
        b = (1 - np.cos(angle)) / square
    return np.eye(3, dtype=dtype) + a * skew + b * (skew @ skew)


def log_rotation(rotation):
    """Branche des petites rotations utilisée par ce diagnostic seulement."""
    dtype = rotation.dtype.type
    w = np.array([rotation[2, 1] - rotation[1, 2],
                  rotation[0, 2] - rotation[2, 0],
                  rotation[1, 0] - rotation[0, 1]], dtype=dtype) / 2
    cosine = (np.trace(rotation) - 1) / 2
    if cosine < -.99:
        raise ValueError('Diagnostic hors domaine : rotation proche de pi')
    sine = np.sqrt(w @ w + dtype(1e-32))
    angle = np.arctan2(sine, cosine)
    factor = (1 + angle * angle / 6 + angle**4 * dtype(7. / 360.)
              if abs(angle) < 1e-3 else angle / sine)
    return factor * w


def rotations_moyennes(rotations, dtype):
    rotations = np.asarray(rotations, dtype=dtype).reshape(-1, 3, 3)
    return np.array([a @ exp_rotation(log_rotation(a.T @ b) / 2)
                     for a, b in zip(rotations[:-1], rotations[1:])])


def inverse3(matrix):
    """Inverse 3×3 sans conversion cachée de longdouble vers float64."""
    cofactor = np.array([np.cross(matrix[1], matrix[2]),
                         np.cross(matrix[2], matrix[0]),
                         np.cross(matrix[0], matrix[1])])
    return cofactor.T / (matrix[0] @ cofactor[0])


def residu_translation(positions, means, lengths, force, cn=CN,
                       deplacements=False, inverse_arrondi=False):
    """Résidu des nœuds libres ; le premier nœud est encastré.

    Le mode déplacements suppose la géométrie de référence rectiligne
    suivant x, avec les mêmes longueurs. Il sépare la déformation de la
    corde de référence AVANT d'y ajouter les petits déplacements.
    inverse_arrondi reproduit le 1/L calculé en f64 par le noyau.
    """
    dtype = means.dtype.type
    positions = np.asarray(positions, dtype=dtype)
    result = np.zeros_like(positions)
    result[-1] = -np.asarray(force, dtype=dtype)
    cn = np.asarray(cn, dtype=dtype)
    for i, (mean, length) in enumerate(zip(means, lengths)):
        inverse_length = (dtype(1. / float(length)) if inverse_arrondi
                          else 1 / dtype(length))
        gamma = (mean.T @ (positions[i + 1] - positions[i])) * inverse_length
        if deplacements:
            reference = mean[0].copy()
            reference[0] -= 1
            gamma += reference
        else:
            gamma[0] -= 1
        internal = mean @ (cn * gamma)
        result[i] -= internal
        result[i + 1] += internal
    return result[1:]


def reconstruit_translations(means, lengths, force, origine, cn=CN):
    """À rotations fixées, chaque poutre porte le même effort F.

    F = R C (Rᵀ d/L - e1) implique
    d = L R⁻ᵀ (e1 + C⁻¹ R⁻¹ F).
    R⁻¹ est utilisé explicitement : les matrices nodales stockées ne sont
    pas remplacées par des rotations réorthonormalisées en haute précision.
    """
    dtype = means.dtype.type
    force = np.asarray(force, dtype=dtype)
    cn = np.asarray(cn, dtype=dtype)
    positions = [np.asarray(origine, dtype=dtype)]
    for mean, length in zip(means, lengths):
        inverse = inverse3(mean)
        stretch = (inverse @ force) / cn
        stretch[0] += 1
        chord = dtype(length) * (inverse.T @ stretch)
        positions.append(positions[-1] + chord)
    return np.array(positions)


def maximum(array):
    return float(np.max(np.abs(array)))


def analyse(ne):
    fraction = .5 * (1 - math.cos(math.pi / 50))
    force = np.array([0., 8.896 * fraction / math.sqrt(2),
                      8.896 * fraction / math.sqrt(2)])
    n, ids = princeton(ne)
    initial = n.etat()
    rest = np.array(initial[1])
    # Même différence f64 que la longueur enregistrée à la construction.
    lengths = np.diff(rest[:, 0])
    n.effort(ids[-1], force.tolist(), [0., 0., 0.])
    usual = appelle(n, False)
    if usual['erreur'] is not None:
        return dict(intervalles=ne, habituel=usual, diagnostic_disponible=False)
    state = n.etat()
    residual = np.array(n.residu_statique()).reshape(-1, 6)[1:]
    means64 = rotations_moyennes(state[2], np.float64)
    means_extended = rotations_moyennes(state[2], np.longdouble)
    same64 = residu_translation(state[1], means64, lengths, force,
                               inverse_arrondi=True)
    same_extended = residu_translation(state[1], means_extended, lengths, force)
    reconstructed = reconstruit_translations(means_extended, lengths, force, rest[0])
    rounded = np.asarray(reconstructed, dtype=np.float64)
    extended_residual = residu_translation(reconstructed, means_extended, lengths, force)
    rounded_residual = residu_translation(rounded, means_extended, lengths, force)
    # Les petits termes doivent être conservés avant l'arrondi en positions
    # absolues : les extraire après cet arrondi ne récupère aucune information.
    displacement = np.asarray(reconstructed - rest.astype(np.longdouble), dtype=np.float64)
    split_residual = residu_translation(displacement, means64, lengths, force,
                                       deplacements=True)
    lost_displacement = rounded - rest
    lost_residual = residu_translation(lost_displacement, means64, lengths, force,
                                      deplacements=True)
    strict_n, strict_ids = princeton(ne)
    strict_n.effort(strict_ids[-1], force.tolist(), [0., 0., 0.])
    strict = appelle(strict_n, True)
    return dict(
        intervalles=ne, diagnostic_disponible=True,
        fraction_charge=fraction, force_bout_N=force.tolist(),
        habituel=usual, strict=strict,
        residus_N=dict(
            noyau_translations=maximum(residual[:, :3]),
            memes_positions_calcul_f64=maximum(same64),
            memes_positions_calcul_etendu=maximum(same_extended),
            reconstruction_etendue=maximum(extended_residual),
            reconstruction_arrondie_f64_calcul_etendu=maximum(rounded_residual),
            deplacements_conserves_f64=maximum(split_residual),
            deplacements_extraits_apres_arrondi=maximum(lost_residual)),
        moments_noyau_Nm=maximum(residual[:, 3:]),
        ecart_forces_noyau_f64_N=maximum(residual[:, :3] - same64),
        ecart_forces_noyau_etendu_N=maximum(residual[:, :3] - same_extended),
        ecart_positions_reconstruction_m=maximum(reconstructed - np.array(state[1], dtype=np.longdouble)),
        pas_representable_x_max_m=maximum(np.spacing(np.asarray(state[1])[:, 0])),
        positions_reference_m=rest.tolist(), longueurs_m=lengths.tolist(),
        etat_habituel=state,
        positions_reconstruites_etendues_m=[[str(x) for x in row] for row in reconstructed],
        deplacements_conserves_f64_m=displacement.tolist(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sortie', type=Path, required=True)
    parser.add_argument('--intervalles', type=int, nargs='+', default=[10, 20, 40, 60])
    args = parser.parse_args()
    if min(args.intervalles) < 1:
        parser.error('Au moins un intervalle requis')
    if np.finfo(np.longdouble).nmant <= np.finfo(np.float64).nmant:
        parser.error('longdouble doit être plus précis que float64 sur cette machine')
    files = ['src/lib.rs', 'src/tangent.rs', 'src/ad.rs', 'ci/confronte_mbdyn.py',
             'ci/diagnostic_statut_statique.py', 'ci/diagnostic_precision_statique.py']
    out = dict(
        perimetre='Premier palier Princeton, formulation milieu ; reconstruction des translations à rotations fixées uniquement',
        python=platform.python_version(), numpy=np.__version__, plateforme=platform.platform(),
        mantisse_bits_hors_bit_implicite={name: np.finfo(dtype).nmant for name, dtype in
                                        [('f64', np.float64), ('longdouble', np.longdouble)]},
        epsilon_longdouble=str(np.finfo(np.longdouble).eps),
        extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
        sources_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        cas=[])
    for ne in args.intervalles:
        case = analyse(ne)
        out['cas'].append(case)
        args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n')
        print(ne, case.get('residus_N', case['habituel']), flush=True)


if __name__ == '__main__':
    main()

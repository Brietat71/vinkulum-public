"""Sépare conditionnement et modèle du câble à rotations bloquées.

Le bloc de translation est assemblé indépendamment avec R C Rᵀ/L.
Aucune source de solveur externe n'est utilisée.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from vinkulum import _vinkulum, campagne
from diagnostic_poutres import modele_factory

ROOT = Path(__file__).resolve().parents[1]


def analyse(formulation):
    parent = modele_factory(formulation)
    class Modele(parent):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.beams = []

        def poutre(self, nom, a, b, ea, ga, gj, ei):
            pa, ra = self.n.pose(a)
            pb, _ = self.n.pose(b)
            chord = np.array(pb)-pa
            length = np.linalg.norm(chord)
            e1 = chord/length
            ya = np.array(ra).reshape(3, 3)[:, 1]
            e2 = ya-e1*np.dot(ya, e1)
            e2 /= np.linalg.norm(e2)
            frame = np.column_stack((e1, e2, np.cross(e1, e2)))
            shear = ga if formulation == 'milieu' else 1/(1/ga+length**2/(12*ei))
            self.beams.append((a, b, frame@np.diag([ea, shear, shear])@frame.T/length))
            return super().poutre(nom, a, b, ea, ga, gj, ei)
    original = campagne.Noyau
    campagne.Noyau = Modele
    try:
        n, ids = campagne._chaine(20, True)
    finally:
        campagne.Noyau = original
    state = n.etat()
    free = np.array([6*i+j for i in ids for j in (0, 2)])
    matrix = np.zeros((6*len(state[1]),)*2)
    for a, b, local in n.beams:
        for i, sign_i in ((a, 1), (b, -1)):
            for j, sign_j in ((a, 1), (b, -1)):
                matrix[6*i:6*i+3, 6*j:6*j+3] += sign_i*sign_j*local
    stiffness = matrix[np.ix_(free, free)]
    kernel = np.array(n.k_c_m_z()[0])[np.ix_(free, free)]
    initial_residual = np.array(n.residu_statique())
    initial_scale = max(1., float(np.max(np.abs(initial_residual[free]))))
    rhs = -initial_residual[free]
    correction = np.linalg.solve(stiffness, rhs)
    positions = np.array(state[1])
    for dof, value in zip(free, correction):
        positions[dof//6, dof % 6] += value
    result = {'tol': 1e-10, 'echelle_forces_initiale': initial_scale,
              'definition_echelle': 'max(1, norme infinie du résidu des translations libres)',
              'conditionnement_2': float(np.linalg.cond(stiffness)),
              'ecart_K_max': float(np.max(np.abs(kernel-stiffness))),
              'ecart_K_relatif': float(np.max(np.abs(kernel-stiffness))/np.max(np.abs(stiffness))),
              'positions_reference': positions.tolist(),
              'fleche_reference_m': float(-positions[:, 2].min()),
              'residu_lineaire_absolu': float(np.max(np.abs(stiffness@correction-rhs))),
              'essais': {}}
    for initialisation in ('originale', 'lineaire'):
        current = list(state)
        if initialisation == 'lineaire':
            current[1] = positions.tolist()
        n.pose_etat(*current)
        start_residual = np.array(n.residu_statique())
        start_scale = max(1., float(np.max(np.abs(start_residual[free]))))
        total_scale = max(1., float(np.max(np.abs(start_residual))))
        start = time.perf_counter()
        try:
            residual, iterations = n.statique(iters=8, paliers_max=1)
            status = {'statut': 'ok', 'residu_annonce': residual, 'iterations': iterations}
        except Exception as e:
            status = {'statut': 'echec', 'erreur': str(e)}
        free_residual = float(np.max(np.abs(np.array(n.residu_statique())[free])))
        status.update(secondes=time.perf_counter()-start,
                      echelle_forces_depart=start_scale,
                      norme_forces_totales_depart=total_scale,
                      residu_relatif_echelle_originale=free_residual/initial_scale,
                      tol_originale_tenue=free_residual/initial_scale <= 1e-10,
                      residu_translations_libres=float(np.max(np.abs(np.array(n.residu_statique())[free]))),
                      phi_max=float(np.max(np.abs(n.phi()))),
                      fleche_m=float(-np.min(np.array(n.etat()[1])[:, 2])))
        result['essais'][initialisation] = status
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sortie', type=Path, required=True)
    args = parser.parse_args()
    result = {'extension_sha256': hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
              'sources_sha256': {str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                                for p in map(Path, ('src/lib.rs','src/tangent.rs',
                                                    'ci/diagnostic_cable.py','ci/diagnostic_poutres.py',
                                                    'python/vinkulum/campagne.py'))},
              'cas': {f:analyse(f) for f in ('milieu','integree')}}
    args.sortie.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(args.sortie)


if __name__ == '__main__':
    main()

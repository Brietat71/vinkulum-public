"""Contrôle de refus, budgets et absence de cache après mutations."""
import sys
import json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
sys.path.insert(0, '/tmp/vinkulum-energie-native-travail/ci')
from trace_complement_dirigee import certifier_complement

r = np.eye(4)
e = np.ones(4)
m = np.eye(4)
b = np.eye(4)[:, :1]
observations = []


def refus(nom, args, kw=None):
    try:
        certifier_complement(*args, **(kw or {}))
    except Exception as ex:
        observations.append(dict(nom=nom, statut='refus', type=type(ex).__name__, message=str(ex)))
    else:
        raise AssertionError(nom+' accepté')


mal = csr_matrix((np.ones(5), np.array([0, 0, 1, 2, 3]), np.array([0, 2, 3, 4, 5])), shape=(4, 4))
mal.has_canonical_format = True
mal.has_sorted_indices = True
refus('D_CSR_doublons_drapeau_mensonger', (mal, r, e, m, b))
refus('R_CSR_doublons_drapeau_mensonger', (r, mal, e, m, b))
refus('M_CSR_doublons_drapeau_mensonger', (r, r, e, mal, b))
lo = r.copy()
lo[1, 0] = .25
refus('R_triangle_inferieur_non_nul', (r, csr_matrix(lo), e, m, b))
ma = m.copy()
ma[1, 0] = .25
refus('M_asymetrique', (r, r, e, csr_matrix(ma), b))
base = certifier_complement(r, r, e, m, b)
rr = r.copy()
rr[0, 0] = 10.
mod = certifier_complement(r, rr, e, m, b)
assert mod['certificat_total']['eta_superieur'] > .98
assert mod['lambda_min'] < .011*base['lambda_min']
observations.append(dict(nom='R_mute_apres_certificat', statut='eta_et_minorant_actualises',
    eta=mod['certificat_total']['eta_superieur'], lambda_min=mod['lambda_min']))
for key in ('budget_operations', 'budget_coefficients', 'budget_rectangulaire', 'precision'):
    refus(key+'_booleen', (r, r, e, m, b), {key: True})
refus('budget_operations_seuil_moins_un', (r, r, e, m, b), {'budget_operations': base['operations_decimal']-1})
assert certifier_complement(r, r, e, m, b, budget_operations=base['operations_decimal'])['lambda_min'] == base['lambda_min']
refus('budget_rectangulaire_seuil_moins_un', (r, r, e, m, b), {'budget_rectangulaire': 4})
assert certifier_complement(r, r, e, m, b, budget_rectangulaire=5)['lambda_min'] == base['lambda_min']
b[1, 0] = 1.
change = certifier_complement(r, r, e, m, b)
assert base['contraintes_sha256'] != change['contraintes_sha256']
observations.append(dict(nom='B_mute_apres_certificat', statut='empreinte_changee_sans_cache'))
Path(__file__).with_suffix('.json').write_text(json.dumps(observations, ensure_ascii=False, indent=2)+'\n')
print(json.dumps(observations, ensure_ascii=False, indent=2))

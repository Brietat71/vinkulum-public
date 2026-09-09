"""Dossier exécutable de certification : décisions exactes, couverture explicite.

L'oracle Fraction utilise as_integer_ratio (CPython), indépendamment du
décodage des bits et de num-bigint (Rust). Une campagne réussie ne transforme
pas les obligations ouvertes en preuves. Aucun seuil statistique d'acceptation.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import struct
import subprocess
import sys
import time

import vinkulum
from vinkulum import _vinkulum
from audit_certification_lineaire import campagne as certificats_lineaires
from audit_certification_quotient import campagne as certificats_quotients


RACINE = Path(__file__).resolve().parents[1]
GRAINE = 0xC37A2026
UNITE = F(1, 1 << 2148)
MARGE = F(1, 1 << 46)
OBLIGATIONS_OUVERTES = {
    'geometrie': 'Étendre les enclosures géométriques 0.15.0 aux fonctions mécaniques et lois temporelles hors du domaine admis.',
    'position': 'Étendre existence locale et distance certifiées en 0.15.0 aux autres mécanismes et grandes dimensions ; garantir la convergence native.',
    'rang': 'Étendre les quotients creux structurés 0.16.0 aux perturbations sans structure de rang et justifier mécaniquement les enveloppes.',
    'dynamique': 'Étendre le certificat initial linéaire aux contacts, partitions et systèmes de chaque pas.',
    'trajectoire': 'Borner l’erreur globale de temps, contraintes, énergie et quantité de mouvement.',
    'contacts': 'Garantir les lois unilatérales, le frottement et les transitions de contact.',
    'flexible': 'Propager les arrondis jusqu’aux réponses physiques des corps flexibles.',
    'implementation': 'Preuve formelle du code Rust/C++ et de ses interfaces ; revue indépendante.',
    'usage': 'Définir et valider expérimentalement les domaines d’usage et leurs incertitudes.',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def flottant(bits):
    return struct.unpack('>d', bits.to_bytes(8, 'big'))[0]


def echantillons():
    petit = flottant(1)
    grand = sys.float_info.max
    u = math.nextafter(1., 0.)
    yield 'frontiere_refusee', [1.], [u], 0., u*(1.-float(MARGE))
    yield 'frontiere_egale', [1.], [1.], 0., 1.-float(MARGE)
    yield 'sous_debordement', [petit], [petit], 0., 0.
    yield 'debordements_annules', [grand, -grand], [grand, grand], 0., 0.
    yield 'annulation_longue', [grand]*512+[-grand]*512+[petit], [grand]*1024+[petit], 0., 0.
    yield 'vide', [], [], -0., 0.
    # Tous les champs d'exposant finis, deux mantisses extrêmes, deux signes.
    # Contrôle de conversion via dt ET de multiplication par un exposant varié.
    for exposant in range(2047):
        for fraction in (0, (1 << 52)-1):
            for signe in (0, 1):
                x = flottant((signe << 63) | (exposant << 52) | fraction)
                y = flottant((((exposant*73+19) % 2047) << 52) | 0x123456789abcd)
                yield 'exposants', [x, -x], [y, 1.], x, abs(x)
    rng = random.Random(GRAINE)

    def reel():
        while True:
            x = flottant(rng.getrandbits(64))
            if math.isfinite(x):
                return x

    for k in range(2048):
        n = rng.randrange(1, 33)
        g, v = [reel() for _ in range(n)], [reel() for _ in range(n)]
        if k % 2 == 0:
            # Beaucoup de décisions positives, même si des produits débordent.
            g, v = g + [-x for x in g], v + v
        yield 'vecteurs', g, v, reel(), abs(reel())
    # Frontières : le budget flottant peut arrondir des deux côtés.
    for _ in range(1024):
        x = flottant((1023 << 52) | rng.getrandbits(52))
        tol = x*(1.-float(MARGE))
        for t in (math.nextafter(tol, 0.), tol, math.nextafter(tol, math.inf)):
            yield 'frontieres', [1.], [x], 0., t


def contre_epreuves():
    compte, acceptes, refuses = {}, 0, 0
    empreinte = hashlib.sha256()
    temoins = []
    debut = time.perf_counter()
    for nom, g, u, dt, tol in echantillons():
        decision, r_hex, s_hex, t_hex = _vinkulum._certificat_ligne_vitesse(g, u, dt, tol)
        produits = [F(a)*F(b) for a, b in zip(g, u, strict=True)]
        r = F(dt)+sum(produits, F(0))
        s = abs(F(dt))+sum(map(abs, produits), F(0))
        attendu = abs(r) <= F(tol)+MARGE*s
        if (int(r_hex, 16)*UNITE, int(s_hex, 16)*UNITE, int(t_hex, 16)*UNITE) != (r, s, F(tol)):
            raise AssertionError(f'{nom}: témoin entier différent de la référence rationnelle')
        if decision != attendu:
            raise AssertionError(f'{nom}: décision exacte incorrecte')
        compte[nom] = compte.get(nom, 0)+1
        acceptes += decision
        refuses += not decision
        entree = {'famille': nom, 'g': [x.hex() for x in g], 'u': [x.hex() for x in u],
                  'dt': dt.hex(), 'tol': tol.hex(), 'accepte': decision,
                  'residu_entier_hex': r_hex, 'echelle_entiere_hex': s_hex,
                  'tolerance_entiere_hex': t_hex}
        empreinte.update(json.dumps(entree, sort_keys=True).encode())
        empreinte.update(b'\n')
        if nom not in ('exposants', 'vecteurs', 'frontieres', 'annulation_longue'):
            temoins.append(entree)
    invalides = [([1.], [], 0., 0.), ([math.nan], [0.], 0., 0.),
                 ([0.], [math.inf], 0., 0.), ([], [], -math.inf, 0.),
                 ([], [], 0., math.nan), ([], [], 0., math.inf), ([], [], 0., -1.)]
    for args in invalides:
        try:
            _vinkulum._certificat_ligne_vitesse(*args)
        except ValueError:
            continue
        raise AssertionError('une entrée hors contrat a été acceptée')
    if not (acceptes and refuses):
        raise AssertionError('contre-épreuves privées de décisions des deux signes')
    return {'statut': 'conforme_sur_la_campagne', 'cas_par_famille': compte,
            'acceptes': acceptes, 'refuses': refuses, 'entrees_invalides_refusees': len(invalides),
            'divergences_oracle': 0, 'graine': GRAINE,
            'empreinte_campagne_sha256': empreinte.hexdigest(), 'temoins': temoins,
            'secondes': time.perf_counter()-debut}


def sources():
    fichiers = subprocess.check_output(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=RACINE
    ).decode().split('\0')
    selection = [p for p in fichiers if p and (
        p.startswith(('src/', 'python/', 'ci/')) or
        p in ('Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'build.rs',
              'docs/CERTIFICATION_NOYAU.md', 'docs/CERTIFICATION_LINEAIRE.md',
              'docs/CERTIFICATION_QUOTIENT.md'))]
    return {p: sha(RACINE/p) for p in sorted(set(selection))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    p.add_argument('--exiger-couverture-complete', action='store_true')
    args = p.parse_args()
    bilan = {
        'schema': 'vinkulum-certification-1', 'version': vinkulum.__version__,
        'noyau_entier_certifie': False, 'certification_industrielle': False,
        'portee': 'Garde de vitesse exact ; bornes linéaires ; rang, compatibilité et unicité mécanique des quotients exacts admis.',
        'preuve_formelle_implementation': False,
        'obligations_ouvertes': OBLIGATIONS_OUVERTES,
        'environnement': {'python': sys.version, 'plateforme': platform.platform(),
                          'extension': str(Path(_vinkulum.__file__).resolve()),
                          'extension_sha256': sha(_vinkulum.__file__)},
        'sources_sha256': sources(),
        'contre_epreuves': contre_epreuves(),
        'certificats_lineaires': certificats_lineaires(),
        'certificats_quotients': certificats_quotients(),
    }
    args.sortie.write_text(json.dumps(bilan, ensure_ascii=False, indent=2)+'\n')
    c = bilan['contre_epreuves']
    print(f"Décisions exactes : {c['acceptes']+c['refuses']} cas, 0 divergence ; "
          f"{c['entrees_invalides_refusees']} entrées invalides refusées.")
    print(f"Systèmes linéaires : {bilan['certificats_lineaires']['certificats_verifies']} certificats vérifiés indépendamment.")
    print(f"Quotients exacts : {bilan['certificats_quotients']['certificats_verifies']} certificats vérifiés indépendamment ; faux rang refusé.")
    print(f"Noyau entier certifié : NON. {len(OBLIGATIONS_OUVERTES)} obligations ouvertes. Dossier : {args.sortie}")
    if args.exiger_couverture_complete:
        raise SystemExit(2)


if __name__ == '__main__':
    main()

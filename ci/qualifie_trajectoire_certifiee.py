"""Trace de la roue installée et certificats temporels du banc pendule."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import vinkulum
from vinkulum import Noyau

import archive_trajectoire_certifiee as archive
from prototype_trajectoire_certifiee import pas


def generer(den):
    if den not in (64, 128, 256):
        raise ValueError('grille hors contrat')
    n = Noyau([0., 0., -9.81])
    n.corps('pendule', 1., [1., 0., 0., 0., 1., 0., 0., 0., 1.], [0., 0., -1.])
    n.liaison('pivot', None, 0, bloque_r=[0, 2])
    etat = n.etat_precis()
    s, c = math.sin(.5), math.cos(.5)
    etat[1][0] = [-s, 0., -c]
    etat[2][0] = [c, 0., s, 0., 1., 0., -s, 0., c]
    n.pose_etat(*etat)
    initial = n.etat()
    # Une seule simulation : ne pas réinitialiser l'accélération algorithmique
    # à chaque observation. La trace publique donne la partie haute de r.
    trajectoire = [initial]+n.simule(1., 1./den, rho=.9, tous=1)
    if len(trajectoire) != den+1:
        raise ValueError('nombre de nœuds natifs non conforme')
    h = F(1, den)
    x, k = archive.INITIAL, archive.COEFFICIENT
    maxima = dict.fromkeys(archive.NOMS, F(0))
    doc = {'schema': archive.SCHEMA, 'modele': archive.MODELE, 'bits': 128,
           'ordre': 8, 'pas': str(h), 'initial': archive.ecrire_boite(x),
           'coefficient': archive.ecrire_intervalle(k), 'noeuds': []}
    for j, ligne in enumerate(trajectoire):
        if F(ligne[0]) != j*h:
            raise ValueError('date native non conforme : subdivision ou sortie absente')
        domaine = None
        if j:
            x, domaine = pas(x, h, k)
        trace = {nom: [str(F(v)) for v in ligne[i][0]]
                 for i, nom in enumerate(archive.NOMS, 1)}
        erreurs = archive.ecarts(x, trace)
        for nom in archive.NOMS:
            maxima[nom] = max(maxima[nom], erreurs[nom])
        doc['noeuds'].append({'t': str(j*h),
                             'domaine': None if domaine is None else archive.ecrire_boite(domaine),
                             'boite': archive.ecrire_boite(x), 'trace': trace})
    doc['majorants'] = {nom: str(maxima[nom]) for nom in archive.NOMS}
    archive.verifier_document(doc)
    return doc


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sortie', type=Path, required=True)
    p.add_argument('--roue', type=Path, required=True)
    args = p.parse_args()
    args.sortie.mkdir(parents=True, exist_ok=False)
    rapport = {'version': vinkulum.__version__, 'python': sys.version,
               'paquet': str(Path(vinkulum.__file__).resolve()),
               'roue_sha256': hashlib.sha256(args.roue.read_bytes()).hexdigest(), 'cas': {},
               'modules_sha256': {nom: hashlib.sha256(Path(__file__).with_name(nom).read_bytes()).hexdigest()
                  for nom in ('intervalle_temporel.py', 'prototype_trajectoire_certifiee.py',
                              'archive_trajectoire_certifiee.py', 'qualifie_trajectoire_certifiee.py')}}
    for den in (64, 128, 256):
        doc = generer(den)
        fichier = args.sortie/('pendule-'+str(den)+'.json')
        fichier.write_text(json.dumps(doc, indent=2)+'\n')
        resultat = subprocess.run([sys.executable, '-S', archive.__file__, str(fichier)],
                                  capture_output=True, text=True, check=True)
        preuve = json.loads(resultat.stdout)
        rapport['cas'][str(den)] = {'sha256': hashlib.sha256(fichier.read_bytes()).hexdigest(),
                                    'verification_autonome': preuve}
    (args.sortie/'qualification.json').write_text(json.dumps(rapport, indent=2)+'\n')
    print(json.dumps(rapport, indent=2))


if __name__ == '__main__':
    main()

"""Requalifie les diagnostics conservés ; les grands champs NPY sont exclus."""
from copy import deepcopy
from pathlib import Path
import unittest
from archive_inertie_binaire import inventaire,lire_json,octets,sha
from experience_controle_facteurs import qualifier,qualifier_controle

ARCHIVE=Path(__file__).resolve().parent.parent/'docs/bancs/livraison-contrainte-0.12.0'


def verifier_resultat(result):
    if result['version']!='0.12.0' or result['certification_machine_reponses'] is not False:
        raise ValueError('identité ou statut de preuve invalide')
    if not qualifier(result) or not all(qualifier_controle(result)):
        raise ValueError('précision ou contrôle refusé')
    expected=1 if result['mode']=='simple' else 3
    if len(result['certificats'])!=expected:raise ValueError('certificats incomplets')
    for cert in result['certificats']:
        for name in ('masse_complete','complement','comparaison_interieure'):
            if cert[name]['certification_machine'] is not True:raise ValueError('preuve manquante')
        if cert['masse_complete']['comparaison_active']['alpha']!=.49:
            raise ValueError('comparaison de masse complète différente')
    return True


def verifier_archive():
    manifest=inventaire(ARCHIVE)
    for group in ('qualification','refus-environnement','refus-metriques'):
        for name,digest in lire_json(ARCHIVE,group+'/manifest.json').items():
            full=group+'/'+name
            actual=(manifest['exclus_npy'][full]['source_sha256'] if full in manifest['exclus_npy']
                    else sha(octets(ARCHIVE,full)))
            if actual!=digest:raise ValueError('archive différente des octets d’origine')
    rows=lire_json(ARCHIVE,'qualification/bilan.json')
    if [r['dossier'] for r in rows]!=[f'{mode}-n{n}' for mode in ('simple','assemblage') for n in (32,128,512)]:
        raise ValueError('plan incomplet')
    for row in rows:
        prefix='qualification/'+row['dossier']
        if row['code']!=0 or lire_json(ARCHIVE,prefix+'.terminal.json')['code']!=0:
            raise ValueError('processus non terminé correctement')
        result=lire_json(ARCHIVE,prefix+'/resultat.json');verifier_resultat(result)
        if manifest['exclus_npy'][prefix+'/champs.npy']['source_sha256']!=result['champs_sha256']:
            raise ValueError('identité du champ exclu incohérente')
    return len(rows)


class LivraisonContrainte(unittest.TestCase):
    def test_inventaire_et_six_qualifications(self):self.assertEqual(verifier_archive(),6)

    def test_contre_epreuves_de_qualification(self):
        original=lire_json(ARCHIVE,'qualification/simple-n32/resultat.json')
        for change in ('certificat','champ','borne','statut'):
            result=deepcopy(original)
            if change=='certificat':result['certificats'][0]['masse_complete']['certification_machine']=False
            elif change=='champ':result['juge']['operateurs']['masse'][0]=1.
            elif change=='borne':result['retours_par_frequence'][0]['bornes']['masse']['relatives'][0]=1.
            else:result['certification_machine_reponses']=True
            with self.assertRaises(ValueError):verifier_resultat(result)

if __name__=='__main__':unittest.main()

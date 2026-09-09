"""Exclusion des échecs du classement et détection des altérations."""
import copy
import math
from pathlib import Path
import unittest

from archive_confrontation_modale import verifier_document, verifier_archive
from confronte_modes import SCHEMA, TAILLES, MOTEURS, SEUIL, REPETITIONS, DELAI, MEMOIRE

DOSSIER=Path(__file__).resolve().parents[1]/'docs/bancs/confrontation-modale-0.17.0'


def donnees():
    # Fixture synthétique pour le juge ; pas une observation de performance.
    refs={str(n):dict(frequences=[1.,2.,3.],bande_hz=[.5,3.5]) for n in TAILLES}
    data=dict(schema=SCHEMA,seuil=SEUIL,tailles=list(TAILLES),moteurs=list(MOTEURS),
              repetitions=REPETITIONS,delai_s=DELAI,memoire_octets=MEMOIRE,
              references=refs,identites={'vinkulum':{'version':'test'},'exudyn':{'version':'test'}},essais=[])
    for n in TAILLES:
        for moteur in MOTEURS:
            for repetition in range(REPETITIONS+1):
                r=dict(moteur=moteur,taille=n,repetition=repetition,echauffement=repetition==0,
                       secondes_processus=.1,rss_kio=1024,code_retour=0,decision='qualifie',
                       erreur_relative_frequence=0.,partie_reelle_relative=0.)
                if moteur.startswith('mbdyn'):
                    z=[(1+.01j*2*math.pi*f)/(1-.01j*2*math.pi*f) for f in (1.,2.,3.)]
                    r['resultat']=dict(dCoef=.01,alpha=[[x.real,x.imag,1.] for x in z])
                else:
                    r['resultat']=dict(version='test',valeurs_propres=[(2*math.pi*f)**2 for f in (1.,2.,3.)])
                data['essais'].append(r)
    return data


class Classement(unittest.TestCase):
    def setUp(self):
        self.d=donnees()

    def verifier(self):
        return verifier_document(self.d,self.d['references'])

    def test_grille_synthetique(self):
        self.assertEqual(sum(x['classe'] for x in self.verifier()),35)

    def test_un_seul_echec_exclut_meme_tres_rapide(self):
        r=self.d['essais'][0]
        r.pop('resultat');r.update(decision='echec_execution',code_retour=1,secondes_processus=1e-6)
        result=self.verifier()[0]
        self.assertFalse(result['classe'])
        self.assertNotIn('mediane_s',result)

    def test_faux_succes_et_spectre_altere(self):
        self.d['essais'][0]['resultat']['valeurs_propres'][0]*=1.001
        with self.assertRaisesRegex(ValueError,'décision spectrale'):
            self.verifier()

    def test_manquant_et_duplique(self):
        row=self.d['essais'].pop()
        with self.assertRaisesRegex(ValueError,'incomplète'):
            self.verifier()
        self.d['essais'].append(copy.deepcopy(self.d['essais'][0]))
        with self.assertRaisesRegex(ValueError,'dupliqué'):
            self.verifier()

    def test_faux_echauffement_et_temps_non_fini(self):
        r=self.d['essais'][0]
        r['echauffement']=False
        with self.assertRaises(ValueError):self.verifier()
        r['echauffement']=True;r['secondes_processus']=float('nan')
        with self.assertRaises(ValueError):self.verifier()

    def test_derniers_bits_de_la_reference(self):
        refs=copy.deepcopy(self.d['references'])
        for ref in refs.values():
            ref['frequences']=[f*(1+2e-14) for f in ref['frequences']]
        result=verifier_document(self.d,refs)
        self.assertTrue(all(r['classe'] for r in result))

    def test_archive_livree(self):
        self.assertEqual(len(verifier_archive(DOSSIER)),35)


if __name__=='__main__':
    unittest.main()

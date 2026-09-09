"""Une décision archivée ne suffit pas si les critères ou rattachements changent."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from defauts_moments_em import encoder
from archive_attribution_reperes import judge,verifier,verifier_observations,sha,rejouer_cas

ARCHIVE=Path(__file__).resolve().parents[1]/'docs/bancs/attribution-reperes-0.18.0'


class ArchiveAttribution(unittest.TestCase):
    def fixture(self):
        budget=encoder(F(64,2**52));zero=encoder(F(0))
        sources={'ci/attribue_moments_em.py':'a'*64,'ci/attribue_rotation_em.py':'b'*64}
        traces={'initial-0.00100.json.gz':'c'*64,'spatial-0.00100.json.gz':'d'*64}
        rv=dict(programme_sha256=sources['ci/attribue_rotation_em.py'],pas_verifies=20000,
                borne_erreur_ode_continue=False,determinants_positifs=True,budget_local=budget,
                budget_rotation=encoder(F(1,10**20)),budget_omega=encoder(F(1,10**18)),verifie=True)
        for key in ('defaut_rotation_local_max','defaut_omega_local_max','orthogonalite_native_max',
                    'orthogonalite_reperes_max','largeur_rotation_max','largeur_omega_max'):rv[key]=zero
        d=dict(cadre='spatial',pas_demande=.001,complet=True,pas_total=20000,pas_verifies=20000,
               sources_sha256=list(traces.values()),programme_sha256=sources['ci/attribue_moments_em.py'],
               budget_local=budget,budget_largeur=encoder(F(1,10**20)),precision_bits=160,
               defauts_locaux_max=[zero,zero],largeur_max=zero,rotation_vitesse=rv,
               attribution_moment_verifiee=True,attribution_rotation_vitesse_verifiee=True)
        return d,traces,sources

    def test_decision_coherente(self):
        d,t,s=self.fixture();judge(d,'rotation-spatial-0.00100',t,s)

    def test_decisions_falsifiees_refusees(self):
        d,t,s=self.fixture()
        mutations=[('pas_verifies',19999),('sources_sha256',['x','y']),('programme_sha256','x'),
                   ('budget_local',encoder(F(1))),('defauts_locaux_max',[]),
                   ('defauts_locaux_max',[encoder(F(1)),encoder(F(0))]),
                   ('largeur_max',encoder(F(-1))),('attribution_moment_verifiee',False)]
        for key,value in mutations:
            with self.subTest(field=key),self.assertRaises(ValueError):
                bad=copy.deepcopy(d);bad[key]=value;judge(bad,'rotation-spatial-0.00100',t,s)

    def test_rotation_falsifiee_refusee(self):
        d,t,s=self.fixture()
        for key,value in [('borne_erreur_ode_continue',True),('determinants_positifs',False),
                          ('orthogonalite_native_max',encoder(F(1))),('budget_omega',encoder(F(1))),
                          ('largeur_rotation_max',encoder(F(1))),('programme_sha256','x')]:
            with self.subTest(field=key),self.assertRaises(ValueError):
                bad=copy.deepcopy(d);bad['rotation_vitesse'][key]=value
                judge(bad,'rotation-spatial-0.00100',t,s)


class ArchiveReelle(unittest.TestCase):
    def test_rejeu_partiel_sans_modification_archive(self):
        with tempfile.TemporaryDirectory() as td:
            target=Path(td)/'archive';shutil.copytree(ARCHIVE,target)
            def hashes():return {str(p.relative_to(target)):sha(p) for p in target.rglob('*') if p.is_file()}
            before=hashes()
            result=rejouer_cas(target,'rotation-materiel-0.00100',Path(td)/'resultat.json',limite=1)
            self.assertEqual(result['pas_verifies'],1)
            self.assertFalse(result['complet'])
            self.assertFalse(result['attribution_rotation_vitesse_verifiee'])
            self.assertEqual(before,hashes())
            self.assertEqual(verifier(target)['cas'],9)

    def test_integrite_et_observations(self):
        self.assertEqual(verifier(ARCHIVE)['cas'],9)
        self.assertEqual(verifier_observations(ARCHIVE)['pics_rattaches'],36)

    def test_rattachement_historique(self):
        from compare_historique_reperes import comparer
        historique=ARCHIVE.parent/'fiabilite-0.14.1-preuves.json.gz'
        self.assertTrue(comparer(ARCHIVE/'traces',historique)['toutes_identiques'])

    def test_alterations_meme_avec_empreinte_recalculee(self):
        with tempfile.TemporaryDirectory() as td:
            target=Path(td)/'archive';shutil.copytree(ARCHIVE,target)
            manifest_path=target/'manifest.json'
            manifest=json.loads(manifest_path.read_text())
            relative='resultats/rotation-spatial-0.00100.json'
            result_path=target/relative;original=json.loads(result_path.read_text())
            def save(data):
                result_path.write_text(json.dumps(data))
                manifest['fichiers_sha256'][relative]=sha(result_path)
                manifest_path.write_text(json.dumps(manifest))
            bad=copy.deepcopy(original)
            bad['defauts_locaux_max'][0]=encoder(F(1))
            save(bad)
            with self.assertRaises(ValueError):verifier(target)
            bad=copy.deepcopy(original)
            peak=bad['rotation_vitesse']['pics_tous_pas']['rotation']
            peak['carre']=encoder(F(1))
            save(bad)
            # Le contrôle de décisions seul ne rejoue pas les observations.
            self.assertEqual(verifier(target)['cas'],9)
            with self.assertRaises(ValueError):verifier_observations(target)
            save(original)
            (target/'traces/initial-0.00100.json.gz').unlink()
            with self.assertRaises(ValueError):verifier(target)


if __name__=='__main__':unittest.main()

"""Contrôles du juge : ne pas classer une référence manquante ou instable."""
import unittest
import gzip
import json
from pathlib import Path
import tempfile
import sys
import time
from types import SimpleNamespace
from unittest.mock import patch

from bilan_confrontation import analyse, difference, load_report
from archive_confrontation import compact
from confronte_mbdyn import run as measure


def run(engine, step, value, ok=True):
    return dict(engine=engine,step=step,rep=0,ok=ok,times=[0.,1.],
                values=[[0.],[value]],wall_seconds=1.,rss_kib=100,
                directory='/tmp/test',error='échec volontaire')


class JugeTests(unittest.TestCase):
    def report(self):
        return dict(metadata={},cases={'test':dict(
            spec=dict(steps=[.1],ref=[.01,.005],target=.1),
            runs=[run(e,h,v) for e in ('vinkulum','mbdyn')
                  for h,v in ((.1,1.02),(.01,1.001),(.005,1.))])})

    def test_accord_et_convergence(self):
        r=analyse(self.report())['cases']['test']
        self.assertTrue(r['references_qualified'])
        self.assertIsNotNone(r['vinkulum_best'])

    def test_reference_en_echec_interdit_classement(self):
        report=self.report()
        report['cases']['test']['runs'][-1]['ok']=False
        r=analyse(report)['cases']['test']
        self.assertFalse(r['references_qualified'])
        self.assertIsNone(r['vinkulum_best'])
        self.assertEqual(len(r['failures']),1)

    def test_desaccord_entre_references_interdit_classement(self):
        report=self.report()
        for r in report['cases']['test']['runs']:
            if r['engine']=='mbdyn': r['values'][-1][0]+=1.
        r=analyse(report)['cases']['test']
        self.assertFalse(r['references_qualified'])
        self.assertIsNone(r['mbdyn_best'])

    def test_reference_non_convergee_interdit_classement(self):
        report=self.report()
        report['cases']['test']['runs'][1]['values'][-1][0]=2.
        self.assertFalse(analyse(report)['cases']['test']['references_qualified'])

    def test_grille_decalee_refusee(self):
        a,b=run('a',.1,1.),run('b',.1,1.)
        b['times']=[.1,1.1]
        with self.assertRaises(ValueError): difference(a,b)

    def test_repetitions_manquantes_interdisent_classement(self):
        report=self.report()
        report['metadata']=dict(protocol_version=2,repetitions=3)
        r=analyse(report)['cases']['test']
        self.assertTrue(r['references_qualified'])
        self.assertIsNone(r['vinkulum_best'])
        self.assertFalse(r['configurations'][0]['repetitions_complete'])

    def test_repetitions_completes_autorisent_classement(self):
        report=self.report()
        report['metadata']=dict(protocol_version=2,repetitions=3)
        runs=report['cases']['test']['runs']
        runs.extend(dict(r,rep=rep) for r in list(runs) if r['step']==.1 for rep in (1,2))
        self.assertIsNotNone(analyse(report)['cases']['test']['vinkulum_best'])
        runs[-1]['ok']=False
        self.assertIsNone(analyse(report)['cases']['test']['mbdyn_best'])

    def test_repetition_dupliquee_refusee(self):
        report=self.report()
        report['metadata']=dict(protocol_version=2,repetitions=1)
        report['cases']['test']['runs'].append(dict(report['cases']['test']['runs'][0]))
        with self.assertRaisesRegex(ValueError,'dupliquée'): analyse(report)

    def test_reference_candidate_non_dupliquee(self):
        report=self.report()
        report['cases']['test']['spec']['steps'].append(.01)
        self.assertEqual(len(analyse(report)['cases']['test']['configurations']),6)

    def test_grilles_invalides_refusees(self):
        for times in ([],[0.],[0.,0.],[1.,0.],[0.,float('nan')]):
            a,b=run('a',.1,1.),run('b',.1,1.)
            a['times']=times
            with self.subTest(times=times),self.assertRaises(ValueError): difference(a,b)

    def test_diagnostic_statique_absent_refuse(self):
        report=self.report()
        report['metadata']=dict(protocol_version=2,repetitions=1)
        report['cases']['test']['spec']['formulation']='integree'
        with self.assertRaisesRegex(ValueError,'paliers stricts'): analyse(report)

    def test_references_propres_a_chaque_moteur(self):
        report=self.report()
        case=report['cases']['test']
        case['spec']['ref_by_engine']={'vinkulum':[.01,.005],'mbdyn':[.02,.001]}
        for r in case['runs']:
            if r['engine']=='mbdyn' and r['step'] in (.01,.005):
                r['step']={.01:.02,.005:.001}[r['step']]
        result=analyse(report)['cases']['test']
        self.assertTrue(result['references_qualified'])
        self.assertIsNotNone(result['vinkulum_best'])
        self.assertIsNotNone(result['mbdyn_best'])

    def test_stagnation_statique_refusee_meme_si_ok(self):
        report=self.report()
        report['metadata']=dict(protocol_version=2,repetitions=1)
        case=report['cases']['test']
        case['spec']['formulation']='integree'
        for r in case['runs']:
            if r['engine']=='vinkulum':
                r['static_diagnostics']=[dict(palier_charge=k,statut='tolerance',strict=True,
                    tolerance_finale_atteinte=True,paliers_stagnation=0,tol=1e-8,
                    residu_relatif=1e-10,contraintes=0.) for k in range(1,51)]
        self.assertIsNotNone(analyse(report)['cases']['test']['vinkulum_best'])
        case['runs'][0]['static_diagnostics'][4]['paliers_stagnation']=1
        with self.assertRaisesRegex(ValueError,'paliers stricts'): analyse(report)

    def test_archive_preserve_les_valeurs_et_detecte_alteration(self):
        report=self.report()
        archive=compact(report)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'trajectoires.json.gz'
            path.write_bytes(gzip.compress(json.dumps(archive).encode()))
            self.assertEqual(analyse(report),analyse(load_report(path)))
            next(iter(archive['trajectoires'].values()))['values'][-1][0]+=1.
            path.write_bytes(gzip.compress(json.dumps(archive).encode()))
            with self.assertRaisesRegex(ValueError,'empreinte'): load_report(path)


class ProcessusTests(unittest.TestCase):
    def test_delai_arrete_aussi_le_solveur_enfant_et_conserve_stdout(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            marker=root/'enfant_encore_actif'
            child=f'import time,pathlib; time.sleep(0.8); pathlib.Path({str(marker)!r}).touch()'
            executable=root/'faux_mbdyn'
            executable.write_text(f'#!{sys.executable}\nimport subprocess,sys,time\n'
                f'subprocess.Popen([sys.executable,"-c",{child!r}])\n'
                'print("enfant démarré",flush=True)\ntime.sleep(20)\n')
            executable.chmod(0o755)
            args=SimpleNamespace(sortie=root,mbdyn=executable,benchmarks=root,timeout=.3)
            with patch('confronte_mbdyn.deck'):
                result=measure(args,'six_barres',.01,'mbdyn',0)
            self.assertFalse(result['ok'])
            self.assertTrue(result['timed_out'])
            self.assertIn('enfant démarré',(Path(result['directory'])/'stdout.txt').read_text())
            time.sleep(.9)
            self.assertFalse(marker.exists(),'un solveur abandonné polluerait les mesures suivantes')


if __name__=='__main__': unittest.main()

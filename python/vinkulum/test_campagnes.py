"""Contrôles de l'ordonnancement, indépendants du coût des modèles physiques."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

from ._campagnes import executer, tache

MODULE = 'vinkulum.test_campagnes'


def _sonde(attente=.2, echec=False):
    debut = time.monotonic()
    print('trace de la sonde', flush=True)
    if echec:
        raise RuntimeError('échec volontaire du cas')
    time.sleep(attente)
    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else []
    return dict(debut=debut, fin=time.monotonic(), cpus=cpus,
                rayon=int(os.environ['RAYON_NUM_THREADS']),
                blas=int(os.environ['OPENBLAS_NUM_THREADS']),
                valeur=sum(i*i for i in range(100)))


class Campagnes(unittest.TestCase):
    def execute(self, cas, jobs=2, timeout=20):
        with tempfile.TemporaryDirectory() as d, patch.dict(os.environ, {
            'VINKULUM_BANCS_LOGS': d, 'VINKULUM_BANCS_WORKER': '0',
            'VINKULUM_BANCS_JOBS': str(jobs), 'VINKULUM_BANCS_CPUS': '2',
            'VINKULUM_BANCS_TIMEOUT': str(timeout),
        }), contextlib.redirect_stdout(io.StringIO()):
            erreur = None
            try:
                valeurs = executer(cas, 'sondes')
            except RuntimeError as e:
                erreur, valeurs = str(e), None
            rapport = json.loads(next(Path(d).glob('*/rapport.json')).read_text())
            logs = [Path(c['log']).read_text() for c in rapport['cas']]
            return valeurs, rapport, logs, erreur

    def test_couverture_budget_et_exclusivite(self):
        cas = [tache(n, MODULE, '_sonde', garder=True, exclusif=n == 'exclusif')
               for n in ('a', 'b', 'exclusif')]
        seq, rs, _, erreur = self.execute(cas, jobs=1)
        self.assertIsNone(erreur)
        par, rp, logs, erreur = self.execute(cas)
        self.assertIsNone(erreur)
        self.assertEqual(list(seq), list(par))
        self.assertEqual(rs['prevus'], rp['prevus'])
        self.assertEqual([c['nom'] for c in rp['cas']], rp['prevus'])
        self.assertTrue(all('trace de la sonde' in log for log in logs))
        for n in par:
            self.assertEqual(seq[n]['valeur'], par[n]['valeur'])
            self.assertEqual(par[n]['blas'], 1)
            if par[n]['cpus']:
                self.assertLessEqual(par[n]['rayon'], len(par[n]['cpus']))
        self.assertGreaterEqual(par['exclusif']['debut'], max(par[n]['fin'] for n in ('a', 'b')))
        if rp['jobs'] == 2:
            self.assertLess(max(par[n]['debut'] for n in ('a', 'b')),
                            min(par[n]['fin'] for n in ('a', 'b')))
            self.assertFalse(set(par['a']['cpus']) & set(par['b']['cpus']))
        self.assertLessEqual(rp['cpu_budget'], 2)
        with patch.dict(os.environ, {'RAYON_NUM_THREADS': '1'}):
            borne, _, _, erreur = self.execute(cas[:1], jobs=1)
            self.assertIsNone(erreur)
            self.assertEqual(borne['a']['rayon'], 1)

    def test_echec_et_timeout_font_echouer_la_campagne(self):
        for kw, timeout in (({'echec': True}, 20), ({'attente': 2}, .2)):
            cas = [tache('mauvais', MODULE, '_sonde', **kw)]
            _, rapport, logs, erreur = self.execute(cas, timeout=timeout)
            self.assertIsNotNone(erreur)
            self.assertEqual(rapport['cas'][0]['statut'], 'echec')
            self.assertEqual(rapport['prevus'], ['mauvais'])
            if kw.get('echec'):
                self.assertIn('échec volontaire', logs[0])

    def test_configuration_et_noms_invalides(self):
        cas = [tache('cas', MODULE, '_sonde')]
        with patch.dict(os.environ, {'VINKULUM_BANCS_WORKER': '0'}):
            with self.assertRaises(ValueError):
                executer(cas*2, 'doublons')
            for variable, valeur in (('VINKULUM_BANCS_JOBS', '0'),
                                     ('VINKULUM_BANCS_CPUS', '-1'),
                                     ('VINKULUM_BANCS_TIMEOUT', 'nan')):
                with patch.dict(os.environ, {variable: valeur}), self.assertRaises(ValueError):
                    executer(cas, 'configuration invalide')


if __name__ == '__main__':
    unittest.main()

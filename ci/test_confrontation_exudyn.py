"""Le classement doit refuser une précision ou une convergence non établie."""
from copy import deepcopy
import gzip
import json
from pathlib import Path
import unittest

from confronte_exudyn import analyse

ARCHIVE = Path(__file__).resolve().parents[1]/'docs/bancs/confrontation-exudyn-0.8.2-essais.json.gz'


class ClassementExudyn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(gzip.decompress(ARCHIVE.read_bytes()))

    def setUp(self):
        self.report = deepcopy(self.original)

    def candidate(self):
        return next(r for r in self.report['runs'] if r['variant']=='exudyn_newton' and r['n']==60 and r['role']=='candidate' and r['rep']==0)

    def test_archive_complete_qualifie_les_references(self):
        r=analyse(self.report)
        self.assertTrue(r['references_qualified'])
        self.assertTrue(all(r['best'].values()))

    def test_palier_non_converge_refuse(self):
        self.candidate()['static_diagnostics'][20]['converged']=False
        with self.assertRaises(ValueError):
            analyse(self.report)

    def test_bilan_force_non_fini_refuse(self):
        self.candidate()['force_balance_max']=float('nan')
        with self.assertRaises(ValueError):
            analyse(self.report)

    def test_binaire_fast_doit_etre_charge(self):
        self.candidate()['exudyn_module']='exudyn.exudynCPPfast'
        with self.assertRaises(ValueError):
            analyse(self.report)

    def test_repetition_dupliquee_refusee(self):
        self.report['runs'].append(deepcopy(self.candidate()))
        with self.assertRaises(ValueError):
            analyse(self.report)

    def test_repetition_manquante_exclut_le_reglage(self):
        self.report['runs'].remove(self.candidate())
        r=analyse(self.report)
        c=next(c for c in r['configurations'] if c['variant']=='exudyn_newton' and c['n']==60)
        self.assertFalse(c['eligible'])

    def test_reference_en_desaccord_exclut_le_classement(self):
        r=next(r for r in self.report['runs'] if r['variant']=='mbdyn' and r['n']==160 and r['role']=='reference' and r['rep']==0)
        r['values'][0][0]+=1e-3
        result=analyse(self.report)
        self.assertFalse(result['references_qualified'])
        self.assertTrue(all(b is None for b in result['best'].values()))

    def test_sensibilite_aux_tolerances_excessive_exclut_le_classement(self):
        r=next(r for r in self.report['runs'] if r['role']=='tolerance')
        r['values'][0][0]+=1e-6
        result=analyse(self.report)
        self.assertFalse(result['references_qualified'])
        self.assertTrue(all(b is None for b in result['best'].values()))


if __name__=='__main__':
    unittest.main()

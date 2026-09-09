"""L'inventaire doit conserver les absences, les écarts et les identités."""
import tempfile
from pathlib import Path
import unittest

from qualifie_validation import historique, rapproche


class Inventaire(unittest.TestCase):
    def test_une_reference_absente_ne_disparait_pas_du_bilan(self):
        ancien = [dict(cas='A', source='mesure', grandeur='force', verdict='ÉCART'),
                  dict(cas='B', source='équation', grandeur='angle', verdict='OK')]
        actuel = [dict(cas='A', source='mesure', grandeur='force', statut='ecart')]
        suivi = rapproche(ancien, actuel)
        self.assertEqual(len(suivi), 2)
        self.assertEqual([r['suivi'] for r in suivi], ['ecart', 'absent'])
        self.assertIsNone(suivi[1]['actuel'])

    def test_les_hypotheses_ne_sont_pas_confondues(self):
        ancien = [dict(cas='Maryland', source='mesure', grandeur=g) for g in ('Ncrit 5', 'Ncrit 9')]
        actuel = [dict(**ancien[1], statut='ok'), dict(**ancien[0], statut='ecart')]
        self.assertEqual([r['suivi'] for r in rapproche(ancien, actuel)], ['ecart', 'ok'])
        with self.assertRaises(ValueError):
            rapproche(ancien, actuel + actuel[:1])

    def test_rapport_historique_et_doublon(self):
        ligne = '| A | mesure | force | 1.0000 | 1.2000 | +20 % | ÉCART | — | note |\n'
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'rapport.md'
            p.write_text('# Rapport\n' + ligne)
            self.assertEqual(historique(p)[0]['verdict'], 'ÉCART')
            p.write_text(ligne*2)
            with self.assertRaises(ValueError):
                historique(p)


if __name__ == '__main__':
    unittest.main()

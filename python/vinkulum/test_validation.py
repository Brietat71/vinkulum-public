"""Une référence en échec doit rester identifiable et faire échouer la commande."""
import contextlib
import io
import unittest
from unittest.mock import patch

from vinkulum import andrews, validation as v, verification as refs

REFERENCES = ('kapitza', 'bennett', 'srscm', 'six_barres', 'lateral_buckling',
              'arbre_tournant', 'quatre_barres_flexible')


class Validation(unittest.TestCase):
    def setUp(self):
        self.lignes = list(v.LIGNES)
        v.LIGNES.clear()
        self.addCleanup(lambda: v.LIGNES.__setitem__(slice(None), self.lignes))

    def campagne(self, echec=None, andrews_echec=None, rapide=False):
        with contextlib.ExitStack() as pile:
            mocks = {nom: pile.enter_context(patch.object(refs, nom,
                     side_effect=RuntimeError('référence indisponible') if nom == echec else None))
                     for nom in REFERENCES}
            pile.enter_context(patch.object(andrews, 'demo', side_effect=andrews_echec))
            v.cas_references(rapide)
        return mocks

    def test_exception_andrews_est_une_erreur_execution(self):
        self.campagne(andrews_echec=RuntimeError('calcul interrompu'))
        self.assertEqual(len(v.LIGNES), 8)
        derniere = v.LIGNES[-1]
        self.assertEqual(derniere[1], 'Schiehlen 1990 (benchmark IUTAM)')
        self.assertFalse(derniere[6])
        self.assertTrue(derniere[7].startswith('ERREUR RuntimeError'))
        # Même décision que le point d'entrée de la commande.
        self.assertTrue(any(l[7].startswith('ERREUR') for l in v.LIGNES))

    def test_une_reference_en_echec_ne_supprime_pas_les_suivantes(self):
        mocks = self.campagne(echec='kapitza')
        self.assertTrue(all(m.call_count == 1 for m in mocks.values()))
        self.assertEqual(len(v.LIGNES), 8)
        self.assertIn('Kapitza', v.LIGNES[0][0])
        self.assertFalse(v.LIGNES[0][6])
        self.assertEqual(sum(l[6] for l in v.LIGNES), 7)

    def test_rapide_conserve_son_perimetre(self):
        mocks = self.campagne(rapide=True)
        self.assertEqual(len(v.LIGNES), 6)
        self.assertEqual(mocks['arbre_tournant'].call_count, 0)
        self.assertEqual(mocks['quatre_barres_flexible'].call_count, 0)

    def test_ecart_physique_ne_devient_pas_exception(self):
        with contextlib.redirect_stdout(io.StringIO()):
            v._ligne('mesure', 'source', 'force', 10., 12., .1)
        self.assertFalse(v.LIGNES[0][6])
        self.assertFalse(v.LIGNES[0][7].startswith('ERREUR'))
        self.assertEqual(v.LIGNES[0][5], .2)


if __name__ == '__main__':
    unittest.main()

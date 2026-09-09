"""Références analytiques du diagnostic d'arrondi, indépendantes de Newton."""
import unittest

import numpy as np

from diagnostic_precision_statique import (
    reconstruit_translations, residu_translation, rotations_moyennes,
)


class PrecisionStatiqueTests(unittest.TestCase):
    def test_chaine_anisotrope_solution_affine(self):
        # Géométrie, coefficients, charges et solution sont binaires exacts.
        for dtype in (np.float64, np.longdouble):
            means = np.repeat(np.eye(3, dtype=dtype)[None], 3, axis=0)
            lengths = np.array([.125, .25, .125], dtype=dtype)
            force, cn = [1., .5, -.25], [8., 4., 2.]
            positions = reconstruit_translations(means, lengths, force, [0., 0., 0.], cn)
            x = np.array([0., .125, .375, .5], dtype=dtype)
            expected = x[:, None] * np.array([1.125, .125, -.125], dtype=dtype)
            np.testing.assert_array_equal(positions, expected)
            np.testing.assert_array_equal(residu_translation(positions, means, lengths, force, cn), 0.)

    def test_meme_solution_tournee_et_translatee(self):
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        means = np.repeat(rotation[None], 2, axis=0)
        force, cn = rotation @ [1., .5, -.25], [8., 4., 2.]
        origine = np.array([3., -2., .5])
        positions = reconstruit_translations(means, [.25, .25], force, origine, cn)
        expected = origine + np.array([0., .25, .5])[:, None] * (rotation @ [1.125, .125, -.125])
        np.testing.assert_array_equal(positions, expected)
        np.testing.assert_array_equal(residu_translation(positions, means, [.25, .25], force, cn), 0.)

    @unittest.skipUnless(np.finfo(np.longdouble).nmant > 52, 'Précision étendue requise')
    def test_extension_perdue_par_position_absolue_f64(self):
        # Barre de longueur 1 : u=F/EA=2^-60, conservable en longdouble
        # et comme déplacement f64, perdu dans la position f64 1+u.
        force, cn = [2.**-30, 0., 0.], [2.**30, 1., 1.]
        means = np.eye(3, dtype=np.longdouble)[None]
        positions = reconstruit_translations(means, [1.], force, [0., 0., 0.], cn)
        reference = np.array([[0., 0., 0.], [1., 0., 0.]])
        self.assertEqual(positions[1, 0] - 1, np.longdouble(2.**-60))
        np.testing.assert_array_equal(residu_translation(positions, means, [1.], force, cn), 0.)
        rounded = positions.astype(np.float64)
        np.testing.assert_array_equal(rounded, reference)
        absolute_residual = residu_translation(rounded, means, [1.], force, cn)
        self.assertEqual(abs(absolute_residual[0, 0]), force[0])
        displacement = np.asarray(positions - reference, dtype=np.float64)
        split_residual = residu_translation(displacement, means.astype(np.float64), [1.], force, cn,
                                           deplacements=True)
        np.testing.assert_array_equal(split_residual, 0.)
        lost_residual = residu_translation(rounded - reference, means.astype(np.float64), [1.], force, cn,
                                          deplacements=True)
        self.assertEqual(abs(lost_residual[0, 0]), force[0])

    def test_rotations_moyennes_plane(self):
        def rz(theta):
            c, s = np.cos(theta), np.sin(theta)
            return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
        angles = np.array([0., .00001, .2, .6])
        means = rotations_moyennes([rz(t) for t in angles], np.float64)
        for mean, theta in zip(means, (angles[1:] + angles[:-1]) / 2):
            np.testing.assert_allclose(mean, rz(theta), rtol=0., atol=3e-16)


if __name__ == '__main__':
    unittest.main()

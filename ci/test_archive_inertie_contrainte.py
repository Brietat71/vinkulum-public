"""Intégrité, requalification et falsifications du lot de 84 essais."""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse import csr_matrix

from archive_inertie_contrainte import (
    ARCHIVE, _ecrire, _experience, _json, _manifest, ecrire_manifest,
    verifier_complet, verifier_inertie, verifier_pivots)
from inertie_complement_dirigee import certifier_inertie_complement


class RejeuInertie(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = csr_matrix(np.diag([1., 1024., 2048.]))
        cls.m = csr_matrix(np.eye(3))
        cls.b = np.array([[1.], [.125], [0.]])
        cls.omega = np.array([0., 2*np.pi*40.])
        cls.parametres = _experience().parametres_certificat("inertie_champs")
        cls.cert = certifier_inertie_complement(cls.d, cls.m, cls.b,
            cls.parametres["gamma_binary64"], precision=32)

    def verifier(self, cert=None, cache=None, b=None):
        verifier_inertie(self.cert if cert is None else cert, self.d, self.m,
            self.b if b is None else b, self.omega, self.parametres, {} if cache is None else cache)

    def test_certificat_rejoue_et_contrainte_exacte(self):
        # ker(B.T): x0=-x1/8. Les deux coefficients du quotient sont
        # strictement positifs à gamma ; aucune valeur propre flottante.
        gamma = Fraction(self.parametres["gamma_binary64"])
        self.assertGreater(Fraction(1024**2)+Fraction(1, 64)-gamma*Fraction(65, 64), 0)
        self.assertGreater(Fraction(2048**2)-gamma, 0)
        self.verifier()

    def test_cache_rejoue_seulement_une_identite(self):
        cache = {}
        self.verifier(cache=cache)
        self.assertEqual(len(cache), 1)
        autre = deepcopy(self.cert)
        autre["phases_s"]["assemblage"] += .01
        autre["preparation_s"] += .01
        with patch("inertie_complement_dirigee.certifier_inertie_complement",
                   side_effect=AssertionError("rejeu inutile")):
            self.verifier(autre, cache)

    def test_faux_pivot_de_meme_signe_refuse_par_rejeu(self):
        cert = deepcopy(self.cert)
        p = cert["preuve_kkt"]["pivots"][0]
        self.assertEqual(len(p["indices"]), 1)
        # La signature déclarée reste mathématiquement valide pour ce
        # faux pivot. Seul le rejeu le rattache réellement à D/M/B.
        p["diagonal"] = [str(Decimal(v)*2) for v in p["diagonal"]]
        verifier_pivots(cert["preuve_kkt"], 4, 32)
        with self.assertRaisesRegex(ValueError, "différent du rejeu"):
            self.verifier(cert)

    def test_un_bit_de_b_change_le_probleme(self):
        b = self.b.copy()
        b[0, 0] = np.nextafter(b[0, 0], np.inf)
        with self.assertRaisesRegex(ValueError, "empreinte B/D/M"):
            self.verifier(b=b)

    def test_booleen_ne_remplace_pas_un_rang(self):
        cert = deepcopy(self.cert)
        cert["rang_contraintes"] = True
        with self.assertRaisesRegex(ValueError, "entières requises"):
            self.verifier(cert)

    def test_pivot_double_determinant_dirige(self):
        preuve = dict(signature=[1, 1, 0], largeur_max=0, pic_coefficients=1,
            coefficients_crees=1, pivots=[dict(indices=[0, 1], signature=[1, 1],
                a=["0", "0"], b=["2", "2"], c=["0", "0"], determinant=["-4", "-4"])])
        self.assertEqual(verifier_pivots(preuve, 2, 32), [1, 1, 0])
        preuve["pivots"][0]["determinant"] = ["-5", "-5"]
        with self.assertRaisesRegex(ValueError, "déterminant dirigé"):
            verifier_pivots(preuve, 2, 32)


class ArchiveInertie(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rapport = _json(gzip.decompress((ARCHIVE/"rapport.json.gz").read_bytes()))

    def test_archive_complete_et_rejeu_reel(self):
        audit = {}
        bilan = verifier_complet(audit=audit)
        self.assertEqual(bilan["nombre_essais"], 84)
        self.assertEqual(len(bilan["configurations"]), 21)
        self.assertGreaterEqual(audit["certificats_inertie_distincts_rejoues"], 3)
        self.assertLessEqual(audit["certificats_inertie_distincts_rejoues"], 24)

    def falsifier(self, mutation, motif):
        with TemporaryDirectory(prefix="vinkulum-inertie-falsifiee-") as tmp:
            rep = Path(tmp)/"archive"
            shutil.copytree(ARCHIVE, rep)
            r = deepcopy(self.rapport)
            mutation(r)
            provenance = _json((rep/"provenance-fichiers.json").read_bytes())
            brut = _json((rep/"manifeste-campagne-brute.json").read_bytes())

            def remplacer(nom, objet):
                meta = provenance[nom]
                data = (json.dumps(objet, ensure_ascii=False, indent=2, allow_nan=False)+"\n").encode()
                if meta["encodage"] != "extrait_rapport":
                    (rep/meta["archive"]).write_bytes(gzip.compress(data, mtime=0)
                        if meta["encodage"] == "gzip" else data)
                meta.update(sha256=hashlib.sha256(data).hexdigest(), octets=len(data))
                brut["fichiers_sha256"][nom] = meta["sha256"]
                return meta["sha256"]

            avant = {e["nom"]: e for e in self.rapport["essais"]}
            for e in r["essais"]:
                if e["resultat"] != avant[e["nom"]]["resultat"]:
                    e["resultat_sha256"] = remplacer("essais/"+e["nom"]+"/resultat.json", e["resultat"])
            remplacer("rapport.json", r)
            _ecrire(rep/"provenance-fichiers.json", provenance)
            _ecrire(rep/"manifeste-campagne-brute.json", brut)
            ecrire_manifest(rep)
            _manifest(rep)  # Les seules empreintes acceptent la falsification.
            with self.assertRaisesRegex(ValueError, motif):
                verifier_complet(rep)

    def test_couverture_falsifiee_apres_rehachage(self):
        self.falsifier(lambda r: r["essais"].pop(), "couverture/ordre")

    def test_faux_pivot_coherent_en_signe_apres_rehachage(self):
        def mutation(r):
            p = r["essais"][0]["resultat"]["certificat_complement"]["preuve_kkt"]["pivots"][0]
            self.assertEqual(len(p["indices"]), 1)
            p["diagonal"] = [str(Decimal(v)*2) for v in p["diagonal"]]
        self.falsifier(mutation, "différent du rejeu")

    def test_faux_majorant_relatif_apres_rehachage(self):
        def mutation(r):
            e = next(e["resultat"] for e in r["essais"] if e["resultat"]["variante"] == "inertie_controle")
            e["retours_par_frequence"][0]["bornes"]["masse"]["relatives"][0] = 0.
        self.falsifier(mutation, "majorant relatif incohérent")

    def test_sondes_integrite_et_octets_originaux(self):
        rep = ARCHIVE.parent/"inertie-contrainte-sondes-2026"
        fichiers = _manifest(rep)
        manifest = _json((rep/"manifest.json").read_bytes())
        self.assertEqual(set(fichiers), set(manifest["fichiers"]))
        for nom, meta in manifest["fichiers"].items():
            self.assertIn(meta["compression"], ("gzip", "aucune"))
            if "source_sha256" in meta:
                data = (rep/nom).read_bytes()
                if meta["compression"] == "gzip":
                    data = gzip.decompress(data)
                self.assertEqual(len(data), meta["source_octets"])
                self.assertEqual(hashlib.sha256(data).hexdigest(), meta["source_sha256"])


if __name__ == "__main__":
    unittest.main()

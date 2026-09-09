"""Contre-épreuves du classement et de la provenance, sans chronométrer de solveur."""
from contextlib import redirect_stdout
from copy import deepcopy
import gzip
import io
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest

from mesure_reponses_groupees import analyse, sha, verifier

ARCHIVE = Path(__file__).resolve().parents[1]/"docs/bancs/reponses-groupees-0.11.0"


class ArchiveReponsesGroupees(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(gzip.decompress((ARCHIVE/"essais.json.gz").read_bytes()))

    def setUp(self):
        # Seule la première configuration est altérée ; ne pas dupliquer les
        # millions de valeurs des autres essais à chaque contre-épreuve.
        self.r = dict(self.original)
        self.r["cas"] = list(self.original["cas"])
        self.r["cas"][0] = deepcopy(self.original["cas"][0])
        self.config = self.r["cas"][0]["configurations"][0]
        self.e = self.config["essais"][1]

    def cible(self):
        return analyse(self.r)["configurations"][0]

    def test_archive_reelle_complete(self):
        with redirect_stdout(io.StringIO()):
            verifier(ARCHIVE)
        b = analyse(self.r)
        self.assertEqual(b, json.loads((ARCHIVE/"bilan.json").read_text()))
        self.assertEqual(len(b["configurations"]), 54)
        self.assertTrue(all(c["complet"] for c in b["configurations"]))
        self.assertEqual(sum(len(v["essais"]) for c in self.r["cas"]
                             for v in c["configurations"]), 216)
        precedent = ARCHIVE.parent/"confrontation-ports-exudyn-0.10.0/bilan.json"
        self.assertEqual(sha(precedent), self.r["choix_hcb_sha256"])

    def test_repetition_manquante_exclut_et_dupliquee_refuse(self):
        self.config["essais"].remove(self.e)
        self.assertFalse(self.cible()["eligible"])
        self.assertNotIn("total_s", self.cible())
        self.config["essais"].extend([self.e, self.e])
        with self.assertRaisesRegex(ValueError, "dupliquée"):
            analyse(self.r)

    def test_couverture_cases_et_variantes(self):
        dernier = self.r["cas"].pop()
        with self.assertRaisesRegex(ValueError, "familles"):
            analyse(self.r)
        self.r["cas"].append(dernier)
        self.r["cas"][0]["configurations"].pop()
        with self.assertRaisesRegex(ValueError, "configurations"):
            analyse(self.r)

    def test_mauvaise_version_et_installation(self):
        self.e["moteur"] = "0.11.0"
        with self.assertRaises(ValueError):
            analyse(self.r)
        self.e["moteur"] = "0.10.0"
        self.e["extension_sha256"] = "0"*64
        with self.assertRaisesRegex(ValueError, "installation"):
            analyse(self.r)

    def test_autre_entree_ou_reference(self):
        for cle in ("entree_sha256", "reference_sha256"):
            ancien = self.e[cle]
            self.e[cle] = "0"*64
            with self.assertRaisesRegex(ValueError, "entrées ou référence"):
                analyse(self.r)
            self.e[cle] = ancien

    def test_budget_cpu_et_versions_numeriques(self):
        self.e["cpu"] = [self.r["cpu"]+1]
        with self.assertRaisesRegex(ValueError, "budget"):
            analyse(self.r)
        self.e["cpu"] = [self.r["cpu"]]
        self.e["numpy"] = "autre"
        with self.assertRaisesRegex(ValueError, "versions numériques"):
            analyse(self.r)

    def test_erreur_operateur_et_maximum(self):
        juge = self.e["juge"]
        juge["operateurs"]["masse"][0] = 2e-6
        with self.assertRaisesRegex(ValueError, "maximum"):
            analyse(self.r)
        juge["maxima_operateurs"]["masse"] = 2e-6
        with self.assertRaisesRegex(ValueError, "qualification"):
            analyse(self.r)
        juge["accepte"] = False
        self.assertFalse(self.cible()["eligible"])

    def test_bornes_distinctes_du_juge_et_couverture(self):
        self.e["controle_bornes"]["masse"][0][0] = False
        self.assertTrue(self.cible()["eligible"])
        self.assertFalse(self.cible()["bornes_confrontees"])
        self.e["controle_bornes"]["masse"].pop()
        with self.assertRaisesRegex(ValueError, "couverture des bornes"):
            analyse(self.r)

    def test_compteur_memoire_et_temps(self):
        self.e["pic_image_kib"] = self.e["rss_image_kib"]-1
        with self.assertRaisesRegex(ValueError, "mémoire"):
            analyse(self.r)
        self.e["pic_image_kib"] = self.e["rss_image_kib"]
        self.e["total_s"] += 1.
        with self.assertRaisesRegex(ValueError, "temps total"):
            analyse(self.r)

    def test_oracle_non_qualifie(self):
        self.r["cas"][0]["reference"]["maxima_operateurs"]["masse"] = 1e-5
        with self.assertRaisesRegex(ValueError, "référence non qualifiée"):
            analyse(self.r)

    def test_archive_alteree_et_provenance_meme_si_manifest_recalcule(self):
        with TemporaryDirectory() as tmp:
            dest = Path(tmp)/"archive"
            shutil.copytree(ARCHIVE, dest)
            fichier = dest/"sources/mesure_reponses_groupees.py"
            fichier.write_text(fichier.read_text()+"\n# modification après mesure\n")
            with self.assertRaisesRegex(ValueError, "altérée"):
                verifier(dest)
            manifest = dest/"manifest.json"
            m = json.loads(manifest.read_text())
            m["fichiers_sha256"]["sources/mesure_reponses_groupees.py"] = sha(fichier)
            manifest.write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError, "provenance de source"):
                verifier(dest)
            fichier.unlink()
            del m["fichiers_sha256"]["sources/mesure_reponses_groupees.py"]
            manifest.write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError, "incomplète"):
                verifier(dest)


if __name__ == "__main__":
    unittest.main()

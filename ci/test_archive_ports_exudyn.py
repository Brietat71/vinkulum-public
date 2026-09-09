"""Qualification et intégrité de l'archive réelle de confrontation des ports.

Ces tests lisent les mesures archivées ; ils ne relancent aucun solveur.
Une archive absente est une erreur de livraison, pas un test ignoré.
Les altérations sont faites sur deepcopy ou dans un répertoire temporaire.
"""
from contextlib import redirect_stdout
from copy import deepcopy
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest

from confronte_ports_exudyn import CIBLE, VARIANTES, analyse, verifier


ARCHIVE = (Path(__file__).resolve().parents[1]
           / "docs/bancs/confrontation-ports-exudyn-0.10.0")
SOURCES_REPRODUCTIBLES = {
    "confronte_ports_exudyn.py", "reference_hcb_exudyn.py",
    "oracle_champs_ports.py", "reference_ports_precision.py",
    "experience_reduction_native.py",
}


class ArchivePortsExudyn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(gzip.decompress((ARCHIVE / "essais.json.gz").read_bytes()))

    def setUp(self):
        self.rapport = deepcopy(self.original)

    def configuration(self):
        cas = next(c for c in self.rapport["cas"] if (c["n"], c["max_hz"]) == (32, 2.))
        return next(c for c in cas["configurations"] if c["variante"] == "vinkulum_reduit")

    def mesure(self):
        return next(e for e in self.configuration()["essais"]
                    if (e["role"], e["repetition"]) == ("mesure", 0))

    @staticmethod
    def resume_cible(bilan):
        return next(c for c in bilan["configurations"]
                    if (c["n"], c["max_hz"], c["variante"]) == (32, 2., "vinkulum_reduit"))

    def test_archive_complete_et_bilan_reproductible(self):
        # Le vérificateur recalcule le bilan et les empreintes déclarées.
        # On contrôle également leur présence et leurs liens au rapport :
        # une liste d'empreintes cohérente mais incomplète ne suffit pas.
        with redirect_stdout(io.StringIO()):
            verifier(ARCHIVE)
        manifeste = json.loads((ARCHIVE / "manifest.json").read_text())
        fichiers = manifeste["fichiers_sha256"]
        presents = {p.relative_to(ARCHIVE).as_posix()
                    for p in ARCHIVE.rglob("*") if p.is_file() and p.name != "manifest.json"}
        self.assertEqual(set(fichiers), presents)
        obligatoires = {"essais.json.gz", "bilan.json", "journaux.json.gz"}
        obligatoires.update("sources/" + nom for nom in SOURCES_REPRODUCTIBLES)
        obligatoires.update(f"entrees/n{n}-f{hz}.npz"
                            for n in (32, 128, 512) for hz in (2, 20, 40))
        self.assertTrue(obligatoires <= fichiers.keys(), obligatoires - fichiers.keys())
        self.assertEqual(manifeste["version_vinkulum"], "0.10.0")
        self.assertTrue(SOURCES_REPRODUCTIBLES <= self.rapport["sources_sha256"].keys())
        for nom, empreinte in self.rapport["sources_sha256"].items():
            self.assertEqual(fichiers["sources/" + nom], empreinte, nom)

        journaux = json.loads(gzip.decompress((ARCHIVE / "journaux.json.gz").read_bytes()))
        noms_journaux = set()
        for cas in self.rapport["cas"]:
            nom = f"n{cas['n']}-f{cas['max_hz']:g}"
            self.assertEqual(fichiers[f"entrees/{nom}.npz"], cas["entree_sha256"], nom)
            self.assertEqual(cas["vinkulum"], manifeste["version_vinkulum"])
            self.assertEqual(cas["reference"]["precision"], [70, 90])
            self.assertRegex(cas["reference"]["fichier_sha256"], r"^[0-9a-f]{64}$")
            essais = cas["sondages"] + [e for c in cas["configurations"] for e in c["essais"]]
            for essai in essais:
                journal = (f"{nom}-{essai['variante']}-r{essai['modes_demandes']}-"
                           f"{essai['role']}{essai['repetition']}.log")
                self.assertNotIn(journal, noms_journaux)
                noms_journaux.add(journal)
                self.assertIn(journal, journaux)
                empreinte = hashlib.sha256(journaux[journal].encode("utf-8")).hexdigest()
                self.assertEqual(empreinte, essai["journal_sha256"], journal)
                self.assertEqual(essai["entree_sha256"], cas["entree_sha256"], journal)
                self.assertEqual(essai["reference_sha256"], cas["reference"]["fichier_sha256"], journal)
                if essai["variante"].startswith("vinkulum"):
                    self.assertEqual(essai["extension_sha256"], cas["extension_sha256"], journal)
                    self.assertEqual(essai["moteur"], cas["vinkulum"], journal)

        bilan = analyse(self.rapport)
        self.assertEqual(bilan, json.loads((ARCHIVE / "bilan.json").read_text()))
        attendus = {(n, float(hz), variante)
                    for n in (32, 128, 512) for hz in (2, 20, 40) for variante in VARIANTES}
        observes = {(c["n"], c["max_hz"], c["variante"]) for c in bilan["configurations"]}
        self.assertEqual(observes, attendus)
        self.assertEqual(len(bilan["configurations"]), 54)
        self.assertTrue(all(c["complet"] for c in bilan["configurations"]))
        # Deux témoins déterminés avant archivage, sans supposer qu'une
        # méthode réussit tous les cas ou qu'elle gagne le classement.
        for variante in ("vinkulum_reduit", "lu_corrigee"):
            cible = next(c for c in bilan["configurations"]
                         if (c["n"], c["max_hz"], c["variante"]) == (32, 2., variante))
            self.assertTrue(cible["eligible"], variante)

    def test_campagne_initiale_et_contre_epreuve_preservees(self):
        initiale = ARCHIVE.with_name(ARCHIVE.name+"-initiale")
        sonde = ARCHIVE.with_name("contre-epreuve-hcb-0.10.0")
        with redirect_stdout(io.StringIO()):
            verifier(initiale)
        r = json.loads((sonde/"hcb-refinement-n32-f20.json").read_text())
        for cle, fichier in (("sonde", sonde/"hcb_refine_probe.py"),
                             ("pont_hcb", initiale/"sources/reference_hcb_exudyn.py"),
                             ("juge", initiale/"sources/confronte_ports_exudyn.py")):
            self.assertEqual(hashlib.sha256(fichier.read_bytes()).hexdigest(), r["sources_sha256"][cle])
        cas = next(c for c in self.original["cas"] if (c["n"], c["max_hz"]) == (32, 20.))
        self.assertEqual(r["entree_sha256"], cas["entree_sha256"])
        self.assertEqual(r["reference_sha256"], cas["reference"]["fichier_sha256"])
        self.assertEqual((r["modes"], r["charges"], len(r["frequences_hz"])), (185, 6, 257))
        for projection in r["projections"]:
            for methode in projection["methodes"]:
                self.assertEqual(methode["statut"], "termine")
                qualifie = max(*methode["maxima"].values(), *methode["maxima_operateurs"].values()) <= CIBLE
                self.assertEqual(methode["accepte"], qualifie)
                if methode["methode"] == "spectrale":
                    self.assertFalse(qualifie)
                else:
                    self.assertTrue(qualifie)

    def test_repetition_dupliquee_refuse_et_manquante_exclut(self):
        self.configuration()["essais"].append(deepcopy(self.mesure()))
        with self.assertRaisesRegex(ValueError, "répétition dupliquée"):
            analyse(self.rapport)
        for role in ("mesure", "chauffe"):
            with self.subTest(role_manquant=role):
                self.rapport = deepcopy(self.original)
                essais = self.configuration()["essais"]
                essais.remove(next(e for e in essais if e["role"] == role))
                cible = self.resume_cible(analyse(self.rapport))
                self.assertFalse(cible["complet"])
                self.assertFalse(cible["eligible"])
                self.assertNotIn("total_s", cible)

    def test_controles_directs_et_compteur_memoire(self):
        racine = Path(__file__).resolve().parents[1]
        r = json.loads((ARCHIVE.parent/"controle-reponse-hcb-0.10.0.json").read_text())
        self.assertEqual(r["source_sha256"], hashlib.sha256(
            (racine/"ci/controle_reponse_hcb.py").read_bytes()).hexdigest())
        for nom, empreinte in r["sources_sha256"].items():
            self.assertEqual(empreinte, self.original["sources_sha256"][nom])
        brut = json.dumps(self.original, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
        self.assertEqual(hashlib.sha256(brut.encode()).hexdigest(), r["campagne_sha256"])
        ecartes = {(c["n"], c["max_hz"], c["variante"]) for c in analyse(self.original)["configurations"]
                   if c["variante"].startswith("hcb") and not c["eligible"]}
        self.assertEqual(ecartes, {(c["n"], c["max_hz"], c["variante"]) for c in r["configurations"]})
        for c in r["configurations"]:
            self.assertEqual({m["methode"] for m in c["methodes"]}, {"directe", "spectrale_corrigee"})
            for m in c["methodes"]:
                j = m["juge"]
                self.assertFalse(j["accepte"])
                self.assertGreater(max(j["maxima_operateurs"].values()), CIBLE)
                for nom, courbe in j["operateurs"].items():
                    self.assertEqual(len(courbe), len(c["indices"]))
                    self.assertEqual(max(courbe), j["maxima_operateurs"][nom])
        r = json.loads((ARCHIVE.parent/"controle-rss-processus-0.10.0.json").read_text())
        self.assertEqual(r["source_sha256"], hashlib.sha256(
            (racine/"ci/controle_rss_processus.py").read_bytes()).hexdigest())
        petit, grand = r["observations"]
        self.assertGreaterEqual(grand["enfant"]["rusage_kib"]-petit["enfant"]["rusage_kib"],
                                grand["parent_reserve_octets"]//1024*9//10)
        hwm = int(next(s for s in grand["enfant"]["proc"] if s.startswith("VmHWM:")).split()[1])
        self.assertLess(hwm, grand["enfant"]["rusage_kib"]//2)

    def test_cas_ou_variante_absent_refuse(self):
        self.rapport["cas"].pop()
        with self.assertRaisesRegex(ValueError, "familles de cas incomplètes"):
            analyse(self.rapport)
        self.rapport = deepcopy(self.original)
        self.rapport["cas"][0]["configurations"].pop()
        with self.assertRaisesRegex(ValueError, "configurations incomplètes"):
            analyse(self.rapport)

    def test_entree_ou_reference_differente_refuse(self):
        for cle in ("entree_sha256", "reference_sha256"):
            with self.subTest(empreinte=cle):
                self.rapport = deepcopy(self.original)
                self.mesure()[cle] = "0" * 64
                with self.assertRaisesRegex(ValueError, "entrées ou référence différentes"):
                    analyse(self.rapport)

    def test_maximum_operateur_corrompu_refuse(self):
        # Un maximum déclaré ne peut masquer les valeurs par fréquence.
        juge = self.mesure()["juge"]
        juge["maxima_operateurs"]["deformation"] = max(juge["operateurs"]["deformation"]) + CIBLE
        with self.assertRaisesRegex(ValueError, "maximum opérateur incohérent"):
            analyse(self.rapport)

    def test_combinaison_de_charges_hors_tolerance_exclut_reglage_complet(self):
        avant = analyse(self.rapport)
        self.assertTrue(self.resume_cible(avant)["eligible"])
        juge = self.mesure()["juge"]
        self.assertTrue(all(v <= CIBLE for v in juge["maxima"].values()))
        # Les six charges individuelles restent qualifiées. Seule l'erreur
        # d'opérateur révèle une combinaison insuffisamment précise.
        juge["operateurs"]["deformation"][0] = 2 * CIBLE
        juge["maxima_operateurs"]["deformation"] = max(juge["operateurs"]["deformation"])
        juge["accepte"] = False
        apres = analyse(self.rapport)
        cible = self.resume_cible(apres)
        self.assertTrue(cible["complet"])
        self.assertFalse(cible["eligible"])
        self.assertEqual(cible["erreurs_operateur_max"]["deformation"], 2 * CIBLE)
        self.assertIn("total_s", cible)  # conserver les mesures du réglage écarté
        autres = lambda bilan: [c for c in bilan["configurations"]
                                if (c["n"], c["max_hz"], c["variante"])
                                != (32, 2., "vinkulum_reduit")]
        self.assertEqual(autres(avant), autres(apres))

    def test_fils_ou_affinite_differents_refuses(self):
        for cause in ("fils", "cpu"):
            with self.subTest(cause=cause):
                self.rapport = deepcopy(self.original)
                essai = self.mesure()
                if cause == "fils":
                    essai["fils"]["OPENBLAS_NUM_THREADS"] = "2"
                else:
                    essai["cpu"] = [self.rapport["cpu"], self.rapport["cpu"] + 1]
                with self.assertRaisesRegex(ValueError, "budget de calcul différent"):
                    analyse(self.rapport)

    def test_verifier_detecte_alterations_meme_avec_manifeste_actualise(self):
        # La copie protège l'archive livrée, y compris quand une assertion
        # échoue. Aucune référence Decimal n'est recalculée ici.
        for alteration in ("empreinte", "bilan", "source_absente", "source_remplacee",
                           "entree", "journal"):
            with self.subTest(alteration=alteration), TemporaryDirectory() as temporaire:
                copie = Path(temporaire) / "archive"
                shutil.copytree(ARCHIVE, copie)
                manifeste = copie / "manifest.json"
                contenu_manifeste = json.loads(manifeste.read_text())
                empreintes = contenu_manifeste["fichiers_sha256"]
                if alteration in ("empreinte", "bilan"):
                    nom = "bilan.json"
                    if alteration == "empreinte":
                        (copie / nom).write_bytes((copie / nom).read_bytes() + b"\n")
                    else:
                        contenu = json.loads((copie / nom).read_text())
                        contenu["cible"] = 2 * CIBLE
                        (copie / nom).write_text(json.dumps(contenu))
                elif alteration.startswith("source_"):
                    nom = "sources/oracle_champs_ports.py"
                    if alteration == "source_absente":
                        (copie / nom).unlink()
                        del empreintes[nom]
                    else:
                        (copie / nom).write_bytes((copie / nom).read_bytes() + b"\n# alteration\n")
                elif alteration == "entree":
                    nom = "entrees/n32-f2.npz"
                    (copie / nom).write_bytes((copie / nom).read_bytes() + b"\n")
                else:
                    nom = "journaux.json.gz"
                    journaux = json.loads(gzip.decompress((copie / nom).read_bytes()))
                    essai = self.mesure()
                    journal = (f"n32-f2-{essai['variante']}-r{essai['modes_demandes']}-"
                               f"{essai['role']}{essai['repetition']}.log")
                    journaux[journal] += "journal altéré\n"
                    (copie / nom).write_bytes(gzip.compress(json.dumps(journaux).encode(), mtime=0))
                if alteration not in ("empreinte", "source_absente"):
                    empreintes[nom] = hashlib.sha256((copie / nom).read_bytes()).hexdigest()
                manifeste.write_text(json.dumps(contenu_manifeste))
                with self.assertRaises(ValueError):
                    verifier(copie)


if __name__ == "__main__":
    unittest.main()

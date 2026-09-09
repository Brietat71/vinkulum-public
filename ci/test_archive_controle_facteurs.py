"""Structure, couverture, rejeu et falsifications du contrôle par facteurs."""
from copy import deepcopy
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest

from archive_controle_facteurs import (
    ARCHIVE, CANDIDATS, VARIANTES, MAILLAGES, _ecrire, _json, _manifest,
    _inventaire, _sha,
    ecrire_manifest, verifier_complet, verifier_facteurs, comparaisons_admises,
)


class DiagnosticsFacteurs(unittest.TestCase):
    def exemple(self):
        u = 2.**-53
        image = dict(rang_facteur=1, colonnes_image=3, defaut_colonnes=[0., 1e-12, 2e-12],
                     kappa=1.000001, certification_machine=False)
        return dict(n=32, directions=3, reduction=dict(taille_complement_reduit=3),
            controle_facteurs=dict(extension="gram", certification_machine=False,
                gamma_gram=192*u/(1-192*u), gamma_sommes=7*u/(1-7*u),
                reparation=dict(masse=deepcopy(image), deformation=deepcopy(image))))

    def test_diagnostic_coherent(self):
        verifier_facteurs(self.exemple())

    def test_gamma_positif_mais_sans_rapport_aux_dimensions(self):
        e = self.exemple()
        e["controle_facteurs"]["gamma_gram"] *= .5
        with self.assertRaisesRegex(ValueError, "différent des dimensions"):
            verifier_facteurs(e)

    def test_rang_booleen_refuse(self):
        e = self.exemple()
        e["controle_facteurs"]["reparation"]["masse"]["rang_facteur"] = True
        with self.assertRaisesRegex(ValueError, "dimensions entières"):
            verifier_facteurs(e)

    def test_defaut_non_fini_et_certification_abusive_refuses(self):
        e = self.exemple()
        e["controle_facteurs"]["reparation"]["masse"]["defaut_colonnes"][0] = float("nan")
        with self.assertRaisesRegex(ValueError, "défauts"):
            verifier_facteurs(e)
        e = self.exemple()
        e["controle_facteurs"]["certification_machine"] = True
        with self.assertRaisesRegex(ValueError, "formule"):
            verifier_facteurs(e)

    def test_aucun_ratio_si_une_configuration_est_refusee(self):
        lignes = [dict(n=n, variante=v, eligible_champs=True,
                      eligible_controle=True, preparation_s=1., reponses_s=2., total_s=3.)
                  for n in MAILLAGES for v in VARIANTES]
        for e in lignes:
            if e["n"] == 32 and e["variante"] == "facteurs_controle":
                e["eligible_controle"] = False
            if e["variante"] == "hcb_standard":
                e["eligible_champs"] = False
            if e["variante"] == "inertie_controle":
                e["eligible_controle"] = False
        valeurs = comparaisons_admises(dict(configurations=lignes))
        self.assertEqual({(e["n"], e["reference"]) for e in valeurs},
                         {(n, v) for n in (128, 512) for v in ("lu_corrigee", "hcb_energie")})


class ArchiveControleFacteurs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rapport = _json(gzip.decompress((ARCHIVE/"rapport.json.gz").read_bytes()))

    def test_archive_complete_et_rejeu_reel(self):
        audit = {}
        bilan = verifier_complet(audit=audit)
        self.assertEqual(bilan["nombre_essais"], 60)
        self.assertEqual(len(bilan["configurations"]), 15)
        self.assertEqual(audit["certificats_inertie_distincts_rejoues"], 3)
        self.assertIs(audit["champs_candidats_sha256_identiques"], True)
        candidats = [e["resultat"] for e in self.rapport["essais"]
                     if e["resultat"]["variante"] in CANDIDATS]
        self.assertEqual(len(candidats), 24)
        self.assertTrue(all(e["controle_oracle"]["accepte_toutes_frequences"] for e in candidats))

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
            e = r["essais"][0]["resultat"]
            e["retours_par_frequence"][0]["bornes"]["masse"]["relatives"][0] = 0.
        self.falsifier(mutation, "majorant relatif incohérent")

    def test_erreur_reelle_depasse_borne_apres_rehachage(self):
        def mutation(r):
            e = r["essais"][0]["resultat"]
            borne = e["retours_par_frequence"][0]["bornes"]["masse"]["absolues"][0]
            norme = e["controle_oracle"]["normes_reference"]["masse"][0][0]
            e["controle_oracle"]["erreurs_absolues"]["masse"][0][0] = 2*borne+1e-8*norme
        self.falsifier(mutation, "majorant confronté incohérent")

    def test_gamma_divise_apres_rehachage(self):
        def mutation(r):
            e = r["essais"][0]["resultat"]
            e["controle_facteurs"]["gamma_gram"] *= .5
        self.falsifier(mutation, "différent des dimensions")


class SondesControleFacteurs(unittest.TestCase):
    def test_integrite_et_octets_originaux(self):
        rep = ARCHIVE.parent/"controle-facteurs-sondes-2026"
        fichiers = _inventaire(rep)
        manifest = _json((rep/"manifest.json").read_bytes())
        self.assertEqual(manifest["schema"], 1)
        self.assertEqual(set(fichiers), set(manifest["fichiers"]))
        for nom, meta in manifest["fichiers"].items():
            self.assertEqual(_sha(fichiers[nom]), meta["sha256"])
            self.assertEqual(fichiers[nom].stat().st_size, meta["octets"])
            self.assertIn(meta["compression"], ("gzip", "aucune"))
            if "source_sha256" in meta:
                data = (rep/nom).read_bytes()
                if meta["compression"] == "gzip":
                    data = gzip.decompress(data)
                self.assertEqual(len(data), meta["source_octets"])
                self.assertEqual(hashlib.sha256(data).hexdigest(), meta["source_sha256"])


if __name__ == "__main__":
    unittest.main()

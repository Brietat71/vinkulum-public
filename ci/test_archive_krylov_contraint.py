"""Contre-épreuves sémantiques après recalcul des empreintes d'archive."""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest

from archive_krylov_contraint import (ARCHIVE, _ecrire, _json, _manifest,
    ecrire_manifest, verifier_complet, verifier_controle)


class ArchiveKrylovContraint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rapport = _json(gzip.decompress((ARCHIVE/"rapport.json.gz").read_bytes()))
        cls.controle = next(e["resultat"] for e in cls.rapport["essais"]
                            if e["resultat"]["variante"] == "krylov_controle")

    def test_archive_integrale(self):
        bilan = verifier_complet()
        self.assertEqual(bilan["nombre_essais"], 60)
        self.assertEqual(len(bilan["configurations"]), 15)

    def falsifier(self, mutation, motif):
        """L'attaquant recalcule les deux manifestes et les SHA des JSON.

        Une incohérence détectée ensuite est donc sémantique, et ne résulte
        pas simplement de l'empreinte d'un fichier édité sans la mettre à jour.
        """
        with TemporaryDirectory(prefix="vinkulum-fausse-archive-") as tmp:
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
                meta["sha256"] = hashlib.sha256(data).hexdigest()
                meta["octets"] = len(data)
                brut["fichiers_sha256"][nom] = meta["sha256"]
                return meta["sha256"]

            # Seuls les résultats modifiés sont réécrits ; aucune source
            # capturée n'est exécutée et aucun NPY de champ n'est requis.
            avant = {e["nom"]: e for e in self.rapport["essais"]}
            for e in r["essais"]:
                if e["resultat"] != avant[e["nom"]]["resultat"]:
                    e["resultat_sha256"] = remplacer("essais/"+e["nom"]+"/resultat.json", e["resultat"])
            remplacer("rapport.json", r)
            _ecrire(rep/"provenance-fichiers.json", provenance)
            _ecrire(rep/"manifeste-campagne-brute.json", brut)
            ecrire_manifest(rep)
            _manifest(rep)  # La couche d'intégrité seule accepte la mutation.
            with self.assertRaisesRegex(ValueError, motif):
                verifier_complet(rep)

    def test_fausse_couverture_apres_rehachage(self):
        self.falsifier(lambda r: r["essais"].pop(), "couverture/ordre")

    def test_pivot_negatif_apres_rehachage(self):
        def mutation(r):
            r["essais"][0]["resultat"]["certificat_complement"]["pivots_h_inferieurs_decimal"] = ["-1"]
        self.falsifier(mutation, "pivots H")

    def test_faux_majorant_relatif_apres_rehachage(self):
        def mutation(r):
            e = next(e["resultat"] for e in r["essais"] if e["resultat"]["variante"] == "krylov_controle")
            e["retours_par_frequence"][0]["bornes"]["masse"]["relatives"][0] = 0.
        self.falsifier(mutation, "majorant relatif incohérent")

    def test_erreurs_absolues_incompatibles_avec_le_juge(self):
        e = deepcopy(self.controle)
        e["controle_oracle"]["erreurs_absolues"]["masse"][0][0] = 0.
        with self.assertRaisesRegex(ValueError, "erreurs absolues et juge relatif"):
            verifier_controle(e)

    def test_majorant_infirmant_oracle(self):
        e = deepcopy(self.controle)
        e["controle_oracle"]["erreurs_absolues"]["masse"][0][0] = 1.
        with self.assertRaisesRegex(ValueError, "majorant confronté incohérent"):
            verifier_controle(e)

    def test_marge_globale_falsifiee(self):
        e = deepcopy(self.controle)
        e["retours_par_frequence"][0]["marge"] *= .5
        with self.assertRaisesRegex(ValueError, "marge globale incohérente"):
            verifier_controle(e)

    def test_archive_sondes_integrite_et_provenance(self):
        rep = ARCHIVE.parent/"krylov-contraint-sondes-2026"
        _manifest(rep)
        # Cette archive exploratoire a un schéma de provenance différent :
        # les SHA/tailles du manifeste couvrent néanmoins tous ses fichiers.


if __name__ == "__main__":
    unittest.main()

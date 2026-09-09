"""Trace les audits archivés des ports relevés, sans lancer de calcul mécanique.

Usage : python ci/trace_ports_releves.py [archive] [préfixe] [--assemblage]
Les maxima d'erreur et les médianes/minima/maxima de temps excluent les
échauffements. Les refus restent visibles et aucun rapport de vitesse n'est
calculé. Les empreintes et la complétude de l'archive sont contrôlées avant
l'import de matplotlib. L'option --assemblage ajoute une figure séparée.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import textwrap


ROOT = Path(__file__).resolve().parents[1]
NOMS = {"direct": "Direct", "qr": "Relèvement QR", "lu_energie": "LU + énergie"}
COULEURS = {"direct": "#916437", "qr": "#087d91", "lu_energie": "#7656a6"}
MARQUEURS = {"direct": "o", "qr": "s", "lu_energie": "D"}
REFUS = "#b13a32"
DUREES = ("preparation_s", "requetes_s", "reconstruction_s")


def lire_json(archive, nom, empreinte=None):
    fichier = (archive / nom).resolve()
    if not fichier.is_relative_to(archive.resolve()):
        raise ValueError(f"Fichier extérieur à l'archive : {nom}")
    contenu = fichier.read_bytes()
    if empreinte is not None and hashlib.sha256(contenu).hexdigest() != empreinte:
        raise ValueError(f"Empreinte incorrecte : {nom}")
    return json.loads(contenu)


def uniques(valeurs, nom):
    valeurs = set(valeurs)
    if len(valeurs) != 1:
        raise ValueError(f"Valeurs incohérentes pour {nom} : {valeurs}")
    return next(iter(valeurs))


def non_negatifs(ligne, champs, nom):
    for champ in champs:
        if not math.isfinite(ligne[champ]) or ligne[champ] < 0:
            raise ValueError(f"{champ} invalide : {nom}")


def lire_archive(archive):
    m = lire_json(archive, "manifest.json")
    if m["schema"] != 1 or set(m["variantes"]) != set(NOMS):
        raise ValueError("Schéma ou variantes non pris en charge")
    cas = {c["nom"]: c for c in m["cas"]}
    if len(cas) != len(m["cas"]) or not cas or m["repetitions"] < 2:
        raise ValueError("Cas dupliqués, cas absents ou répétitions insuffisantes")
    for nom, empreinte in m["sources"].items():
        contenu = (archive / "sources" / nom).read_bytes()
        if hashlib.sha256(contenu).hexdigest() != empreinte:
            raise ValueError(f"Empreinte de source incorrecte : {nom}")
    references = {nom: lire_json(archive, "references/" + nom, empreinte)
                  for nom, empreinte in m["references"].items()}
    groupes = {(nom, v): [] for nom in cas for v in m["variantes"]}
    assemblages = {f"assemblage-{n}-{s}": [] for n, s in m["assemblages"]}
    attendus = {(nom, v, r) for nom, v in groupes for r in range(m["repetitions"])}
    attendus |= {(nom, "assemblage_qr", r) for nom in assemblages
                 for r in range(m["repetitions"])}
    vus, mesures = set(), []
    for e in m["essais"]:
        cle = (e["cas"], e["variante"], e["repetition"])
        if cle in vus or cle not in attendus:
            raise ValueError(f"Essai inconnu ou dupliqué : {cle}")
        vus.add(cle)
        d = lire_json(archive, e["fichier"], e["sha256"])
        if (d["variante"] != e["variante"]
                or e["echauffement"] != (e["repetition"] == 0)
                or d["certification_machine"] is not False):
            raise ValueError(f"Identité, échauffement ou certification incohérents : {cle}")
        non_negatifs(d, (*DUREES, "total_s", "borne_uniforme", "directions", "ports"), cle)
        if d["total_s"] <= 0 or not math.isclose(
                d["total_s"], sum(d[k] for k in DUREES), rel_tol=1e-12):
            raise ValueError(f"Durées incohérentes : {cle}")
        if len(d["cpu"]) != 1 or any(v != "1" for v in d["fils"].values()):
            raise ValueError(f"Configuration CPU incohérente : {cle}")
        if d["variante"] == "assemblage_qr":
            if e["cas"] != f'assemblage-{d["n"]}-{d["segments"]}':
                raise ValueError(f"Assemblage incohérent : {cle}")
            non_negatifs(d, ("erreur_port_max", "erreur_relative_port_max",
                             "rapport_erreur_borne_max", "reponses_non_bornees"), cle)
            destination = assemblages[e["cas"]]
        else:
            non_negatifs(d, ("erreur_D", "erreur_physique", "tolerance"), cle)
            nom_reference = e["cas"] + ".json"
            if (d["cas"] != cas[e["cas"]] or d["tolerance"] <= 0
                    or d["oracle_sha256"] != m["references"][nom_reference]
                    or references[nom_reference]["cas"] != d["cas"]):
                raise ValueError(f"Cas ou référence incohérents : {cle}")
            accepte = (d["borne_uniforme"] <= d["tolerance"]
                       and d["erreur_D"] <= d["tolerance"]
                       and d["erreur_physique"] <= d["tolerance"])
            if d["accepte"] != accepte:
                raise ValueError(f"Acceptation incohérente : {cle}")
            destination = groupes[cle[:2]]
        if not e["echauffement"]:
            destination.append(d)
            mesures.append(d)
    if vus != attendus:
        raise ValueError("Archive incomplète : répétitions manquantes")
    repetitions = uniques((len(ds) for ds in (*groupes.values(), *assemblages.values())),
                          "répétitions hors échauffement")
    normales = [d for ds in groupes.values() for d in ds]
    tolerance = uniques((d["tolerance"] for d in normales), "tolérance")
    for nom, ds in assemblages.items():
        for d in ds:
            accepte = (d["borne_uniforme"] <= tolerance and d["reponses_non_bornees"] == 0
                       and d["rapport_erreur_borne_max"] <= 1)
            if d["accepte"] != accepte:
                raise ValueError(f"Acceptation incohérente : {nom}")
    meta = dict(tolerance=tolerance, repetitions=repetitions,
                requetes=uniques((d["nombre_requetes"] for d in normales), "requêtes"),
                audit=uniques((d["nombre_frequences_audit"] for d in normales), "fréquences d'audit"),
                dps=uniques((r["dps"] for r in references.values()), "précision de l'oracle"),
                cpu=uniques((tuple(d["cpu"]) for d in mesures), "affinité CPU"))
    return m, groupes, assemblages, meta


def nombre(x, chiffres=3):
    return f"{x:.{chiffres}g}".replace(".", ",")


def entier(x):
    return f"{x:,.0f}".replace(",", "\u202f")


def statistiques(ds, champ, facteur=1.):
    valeurs = [d[champ] * facteur for d in ds]
    mediane = statistics.median(valeurs)
    return mediane, mediane - min(valeurs), max(valeurs) - mediane


def plancher(valeurs):
    positifs = [v for v in valeurs if v > 0]
    return min(positifs) / 4 if positifs else 1e-16


def point(ax, position, valeur, accepte, couleur, marqueur, minimum, erreurs=None):
    couleur = couleur if accepte else REFUS
    marqueur = marqueur if accepte else "x"
    affichee = max(valeur, minimum)
    ax.errorbar(position, affichee, yerr=erreurs, fmt=marqueur,
                color=couleur, ecolor=couleur, markersize=7, markeredgewidth=1.6,
                capsize=4, linewidth=1.25, zorder=4)
    if valeur == 0:
        ax.annotate("0", (position, affichee), xytext=(0, -12),
                    textcoords="offset points", ha="center", fontsize=8)


def decorer(ax, labels, titre, unite):
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.set_yscale("log")
    ax.set_xticks(range(len(labels)), labels, fontsize=10)
    ax.tick_params(axis="x", length=0, pad=10)
    ax.set_title(titre, loc="left", fontsize=13, pad=16)
    ax.set_ylabel(unite)
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: nombre(v)))
    ax.grid(axis="y", which="major", color="#dbe2e5", linewidth=.8)
    ax.set_axisbelow(True)
    ax.margins(x=.12, y=.3)


def pied(fig, archive, m, meta, notes):
    try:
        source = archive.resolve().relative_to(ROOT)
    except ValueError:
        source = archive.resolve()
    fig.text(.065, .225, "\n".join(textwrap.wrap(notes, 165)), fontsize=10,
             color="#596970", va="top", linespacing=1.5)
    fig.text(.065, .09,
             f'CPU : {m["cpu_info"]} · cœur logique {meta["cpu"][0]} · un fil par bibliothèque.',
             fontsize=9.5)
    fig.text(.065, .06,
             "Croix rouges : au moins un essai refusé. Aucun rapport de vitesse avec une variante refusée. "
             "Bornes sans certification des arrondis.", fontsize=9.5)
    fig.text(.065, .03, f"Source : {source}/manifest.json et essais bruts vérifiés par SHA256.",
             fontsize=9, color="#687b84")


def tracer_ports(archive, m, groupes, meta):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 8.3))
    fig.subplots_adjust(left=.065, right=.975, bottom=.31, top=.72, wspace=.23)
    fig.suptitle("Ports relevés : précision physique et coût", x=.065, y=.967,
                 ha="left", fontsize=20, fontweight="bold")
    fig.text(.065, .914,
             f'Tolérance du Schur normalisé : {nombre(meta["tolerance"])} · '
             f'{meta["repetitions"]} essais hors échauffement par cas et variante.', fontsize=11)
    fig.text(.065, .873,
             f'Maximum d’erreur sur {meta["audit"]} fréquences auditées × {meta["repetitions"]} essais ; '
             "temps médian et intervalle minimum–maximum.", fontsize=10.5)
    labels = []
    toutes = [d for ds in groupes.values() for d in ds]
    minimum = plancher([meta["tolerance"], *[d["erreur_physique"] for d in toutes]])
    for i, c in enumerate(m["cas"]):
        famille = {"chaine": "Chaîne", "console": "Console"}.get(c["famille"], c["famille"])
        labels.append(f'{famille}\n{entier(c["n"])} ' +
                      ("masses" if c["famille"] == "chaine" else "éléments"))
        for j, variante in enumerate(m["variantes"]):
            ds = groupes[(c["nom"], variante)]
            position = i + (j - (len(m["variantes"]) - 1) / 2) * .23
            accepte = all(d["accepte"] for d in ds)
            point(axes[0], position, max(d["erreur_physique"] for d in ds), accepte,
                  COULEURS[variante], MARQUEURS[variante], minimum)
            mediane, bas, haut = statistiques(ds, "total_s", 1000.)
            point(axes[1], position, mediane, accepte, COULEURS[variante], MARQUEURS[variante],
                  0., [[bas], [haut]])
    axes[0].axhline(meta["tolerance"], color="#555f64", linestyle="--", linewidth=1.2)
    decorer(axes[0], labels, "Erreur physique maximale observée", "Norme 2 du défaut de Schur normalisé")
    decorer(axes[1], labels, "Temps total par problème", "Temps écoulé sur un cœur (ms)")
    handles = [Line2D([], [], color=COULEURS[v], marker=MARQUEURS[v], linestyle="none",
                      markersize=7, label=NOMS[v]) for v in m["variantes"]]
    handles += [Line2D([], [], color=REFUS, marker="x", linestyle="none", markersize=7,
                       label="Refus"),
                Line2D([], [], color="#555f64", linestyle="--", label="Tolérance")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.06, .83),
               ncol=len(handles), frameon=False, fontsize=10.5)
    ordre = ", ".join(NOMS[v] for v in m["variantes"])
    notes = (f'Dans chaque cas, de gauche à droite : {ordre}. '
             f'Coût : construction + contrôle uniforme + {meta["requetes"]} requêtes + un champ. '
             f'Oracles à {meta["dps"]} chiffres, imports et modèles exclus des durées. '
             'Les points d’erreur sont des audits échantillonnés, pas un supremum certifié sur la bande. '
             'Les zéros éventuels sont placés au pied du logarithme et annotés « 0 ».')
    pied(fig, archive, m, meta, notes)
    return fig


def tracer_assemblages(archive, m, groupes, meta):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    if not groupes:
        raise ValueError("Aucun assemblage dans l'archive")
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 8.3))
    fig.subplots_adjust(left=.06, right=.98, bottom=.31, top=.73, wspace=.32)
    fig.suptitle("Ports relevés : réponse aux ports après assemblage", x=.065, y=.967,
                 ha="left", fontsize=20, fontweight="bold")
    requetes = uniques((d["nombre_frequences"] for ds in groupes.values() for d in ds),
                      "requêtes de l'assemblage")
    fig.text(.065, .913,
             f'{meta["repetitions"]} essais hors échauffement · {requetes} réponses par assemblage · '
             f'cible du contrôle de Schur assemblé : {nombre(meta["tolerance"])}.', fontsize=11)
    fig.text(.065, .87,
             "Maxima pour l’erreur ; médianes et minimum–maximum pour les temps et les dimensions. "
             "Annotations εrel : erreur relative maximale aux ports.", fontsize=10.5)
    labels = []
    minimum = plancher([1., *[d["rapport_erreur_borne_max"]
                             for ds in groupes.values() for d in ds]])
    for i, (n, segments) in enumerate(m["assemblages"]):
        ds = groupes[f"assemblage-{n}-{segments}"]
        labels.append(f'{entier(n)} masses\n{entier(segments)} sous-structures')
        accepte = all(d["accepte"] for d in ds)
        rapport = max(d["rapport_erreur_borne_max"] for d in ds)
        point(axes[0], i, rapport, accepte, COULEURS["qr"], "s", minimum)
        erreur_relative = max(d["erreur_relative_port_max"] for d in ds)
        annotation = f"εrel = {nombre(erreur_relative, 2)}"
        sans_borne = max(d["reponses_non_bornees"] for d in ds)
        if sans_borne:
            annotation += f"\n{sans_borne} réponses sans borne"
        axes[0].annotate(annotation, (i, max(rapport, minimum)), xytext=(0, 12),
                         textcoords="offset points", ha="center", fontsize=9)
        mediane, bas, haut = statistiques(ds, "total_s", 1000.)
        point(axes[1], i, mediane, accepte, COULEURS["qr"], "s", 0., [[bas], [haut]])
        for champ, decalage, couleur, marqueur in (
                ("directions", -.12, COULEURS["qr"], "s"),
                ("ports", .12, "#916437", "o")):
            mediane, bas, haut = statistiques(ds, champ)
            if mediane - bas <= 0:
                raise ValueError("Dimension nulle incompatible avec le panneau logarithmique")
            point(axes[2], i + decalage, mediane, accepte, couleur, marqueur, 0., [[bas], [haut]])
            axes[2].annotate(entier(mediane), (i + decalage, mediane + haut), xytext=(0, 7),
                             textcoords="offset points", ha="center", fontsize=9)
    axes[0].axhline(1., linestyle="--", color="#555f64", linewidth=1.2)
    decorer(axes[0], labels, "Réponse aux ports : erreur / borne", "Maximum du rapport en norme 2")
    decorer(axes[1], labels, "Coût de l'assemblage", "Temps écoulé sur un cœur (ms)")
    decorer(axes[2], labels, "Dimensions conservées", "Nombre de directions ou de ports")
    fig.legend(handles=[
        Line2D([], [], color=COULEURS["qr"], marker="s", linestyle="none", label="QR / directions intérieures"),
        Line2D([], [], color="#916437", marker="o", linestyle="none", label="Ports globaux"),
        Line2D([], [], color=REFUS, marker="x", linestyle="none", label="Refus"),
        Line2D([], [], color="#555f64", linestyle="--", label="Erreur / borne = 1")],
        loc="upper left", bbox_to_anchor=(.06, .83), ncol=4, frameon=False, fontsize=10.5)
    notes = (f'Coût : construction des sous-structures et contrôle + {requetes} réponses + un champ. '
             'Référence analytique hors durées. Le rapport porte sur l’erreur absolue globale aux ports '
             'divisée par sa borne, avec maximum sur les fréquences disposant d’une borne ; toute réponse '
             'sans borne entraîne un refus. Les zéros éventuels sont annotés « 0 ».')
    pied(fig, archive, m, meta, notes)
    return fig


def sauver(fig, prefixe, archive):
    import matplotlib.pyplot as plt
    prefixe.parent.mkdir(parents=True, exist_ok=True)
    empreinte = hashlib.sha256((archive / "manifest.json").read_bytes()).hexdigest()
    description = f"Archive {archive}, manifeste SHA256 {empreinte}. Échauffements exclus, refus conservés."
    svg, png = prefixe.with_suffix(".svg"), prefixe.with_suffix(".png")
    fig.savefig(svg, metadata={"Date": None, "Title": "Vinkulum : ports relevés",
                              "Description": description, "Creator": "Vinkulum"})
    svg.write_text("\n".join(s.rstrip() for s in svg.read_text().splitlines()) + "\n")
    fig.savefig(png, dpi=180, metadata={"Description": description, "Software": "Vinkulum"})
    plt.close(fig)
    print(svg)
    print(png)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path,
                        default=ROOT / "docs/bancs/ports-releves-2026")
    parser.add_argument("sortie", nargs="?", type=Path,
                        default=ROOT / "docs/figures/ports-releves-2026",
                        help="préfixe des fichiers .svg et .png")
    parser.add_argument("--assemblage", action="store_true",
                        help="ajoute les fichiers <préfixe>-assemblage.svg/.png")
    args = parser.parse_args()
    try:
        m, groupes, assemblages, meta = lire_archive(args.archive)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({
            "font.family": "DejaVu Sans", "font.size": 11,
            "axes.spines.top": False, "axes.spines.right": False,
            "axes.titleweight": "bold", "axes.labelcolor": "#27343b",
            "text.color": "#27343b", "axes.edgecolor": "#bcc7cc",
            "svg.hashsalt": "vinkulum-ports-releves-2026", "svg.fonttype": "none",
        })
        prefixe = args.sortie.with_suffix("") if args.sortie.suffix in (".svg", ".png") else args.sortie
        figures = [(tracer_ports(args.archive, m, groupes, meta), prefixe)]
        if args.assemblage:
            figures.append((tracer_assemblages(args.archive, m, assemblages, meta),
                            prefixe.with_name(prefixe.name + "-assemblage")))
        for fig, sortie in figures:
            sauver(fig, sortie, args.archive)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()

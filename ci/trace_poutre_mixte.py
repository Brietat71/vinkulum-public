"""Trace les résultats archivés du prototype ; ne relance aucun calcul."""

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefixe", type=Path)
    args = parser.parse_args()
    prefix = str(args.prefixe)
    archive = json.loads(gzip.decompress(Path(prefix + "-donnees.json.gz").read_bytes()))
    campaign = json.loads(archive["campagne_json"])
    diagnostics = json.loads(archive["diagnostics_json"])
    probes = [json.loads(line) for line in diagnostics["linearisation"]["stdout"].splitlines()]
    colors = {1: "#246b91", 2: "#b14e2b"}
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for order in (1, 2):
        rows = [r for r in campaign["configurations"] if r["ordre"] == order and r["reussites"] == 3]
        n = [r["elements"] for r in rows]
        for field, style, label in (("erreur_ligne_m", "-o", "ligne, 961 points"), ("erreur_bout_m", "--s", "bout")):
            axes[0, 0].loglog(n, [r[field] * 1e6 for r in rows], style, color=colors[order], label=f"Ordre {order} · {label}")
        axes[0, 1].loglog([r["temps_median_s"] for r in rows], [r["erreur_ligne_m"] * 1e6 for r in rows], "-o", color=colors[order], label=f"Ordre {order}")
        for r in rows:
            axes[0, 1].annotate(str(r["elements"]), (r["temps_median_s"], r["erreur_ligne_m"] * 1e6), xytext=(4, 5), textcoords="offset points", fontsize=8)
        moment = [r for r in diagnostics["moment"] if r["ordre"] == order and r["code"] == 0]
        axes[1, 0].loglog([r["elements"] for r in moment], [r["erreur_ligne_m"] for r in moment], "-o", color=colors[order], label=f"Ordre {order}")
        linear = [r for r in probes if r["ordre"] == order and r["rigidite_flexion"] > 0]
        axes[1, 1].loglog([r["rigidite_cisaillement"] for r in linear], [max(r["erreur_flexibilite_equilibree"], 1e-16) for r in linear], "-o", color=colors[order], label=f"Ordre {order}")
    for ax in axes[0]:
        ax.axhline(10, color="#555555", linestyle=":", linewidth=1)
    axes[0, 0].set(title="Princeton : le bout ne suffit pas", xlabel="Nombre d'éléments", ylabel="Erreur maximale par composante (µm)")
    axes[0, 1].set(title="Coût du prototype seul · 3 répétitions", xlabel="Temps total médian (s) · numéro = éléments", ylabel="Erreur sur 961 points (µm)")
    axes[1, 0].set(title="Moment pur : comparaison à l'arc exact", xlabel="Nombre d'éléments", ylabel="Erreur sur 961 points (m)")
    axes[1, 1].set(title="Flexibilité linéaire : limite de la condensation", xlabel="GA L² / EI", ylabel="Erreur relative équilibrée")
    axes[1, 1].annotate("Raideur négative\nartificielle", xy=(1e16, 5.0887), xytext=(1e7, 1e-3), arrowprops={"arrowstyle": "->"}, fontsize=9)
    for ax, ticks in ((axes[0, 0], [1, 2, 4, 8, 16]), (axes[0, 1], [0.3, 0.5, 1, 2, 4]), (axes[1, 0], [1, 2, 4, 8])):
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    for ax in axes.flat:
        ax.grid(True, which="major", alpha=0.22)
        ax.legend(fontsize=8)
    fig.suptitle("Vinkulum · prototype statique mixte issu d'une formulation de 2026", fontsize=14)
    for extension in ("svg", "png"):
        path = Path(prefix + "." + extension)
        fig.savefig(path, dpi=160)
        if extension == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")


if __name__ == "__main__":
    main()

"""Figure autonome des coûts à précision physique commune, depuis le bilan."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter


STYLES = {
    "vinkulum_api": ("Vinkulum · API avec bornes", "#15636c", "s"),
    "vinkulum_reduit": ("Vinkulum · 6 charges groupées", "#00a58a", "o"),
    "lu": ("LU creuse équilibrée", "#747b86", "^"),
    "lu_corrigee": ("LU + 2 corrections par D", "#252c36", "v"),
    "hcb_standard": ("HCB Exudyn · projection K", "#d58b22", "D"),
    "hcb_energie": ("HCB Exudyn · projection D", "#a44e8c", "P"),
}


def tracer(bilan, destination):
    donnees = json.loads(Path(bilan).read_text())["configurations"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 5.2), sharey=True)
    handles = []
    for ax, hz in zip(axes, (2, 20, 40)):
        for variante, (label, color, marker) in STYLES.items():
            points = sorted((v for v in donnees if v["max_hz"] == hz and v["variante"] == variante),
                            key=lambda v: v["n"])
            bons = [v for v in points if v["eligible"]]
            ax.plot([v["n"] for v in bons], [1000*v["total_s"] for v in bons],
                    marker=marker, color=color, linewidth=1.4, markersize=6, label=label)
            mauvais = [v for v in points if not v["eligible"] and "total_s" in v]
            ax.scatter([v["n"] for v in mauvais], [1000*v["total_s"] for v in mauvais],
                       marker="x", color=color, s=50, linewidth=1.5)
            if hz == 2:
                handles.append(Line2D([], [], color=color, marker=marker, label=label))
        ax.set_xscale("log", base=2)
        ax.set_xticks([32, 128, 512], ["32", "128", "512"])
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.grid(axis="y", alpha=.22)
        ax.set_title(f"0–{hz} Hz", loc="left", weight="bold")
        ax.set_xlabel("Nombre de poutres")
        if hz == 40:
            ax.text(.03, .97, "Vinkulum : bande refusée", transform=ax.transAxes,
                    va="top", fontsize=9, color="#15636c")
    axes[0].set_ylabel("Préparation + 257 fréquences × 6 charges (ms)")
    handles.append(Line2D([], [], color="#555555", marker="x", linestyle="none",
                          label="Hors seuil de précision : exclu du classement"))
    fig.suptitle("Réduction par ports : coût total et qualification physique", x=.065,
                 ha="left", fontsize=15, weight="bold")
    fig.text(.065, .90, "Médianes de 3 processus frais · un cœur · seuil relatif 10⁻⁶ · champs, déformations et ports",
             fontsize=10, color="#525967")
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(.52, -.02), fontsize=9)
    fig.subplots_adjust(top=.79, bottom=.27, left=.075, right=.98, wspace=.15)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(destination)+".svg", bbox_inches="tight")
    svg = Path(str(destination)+".svg")
    svg.write_text("\n".join(ligne.rstrip() for ligne in svg.read_text().splitlines())+"\n")
    fig.savefig(str(destination)+".png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("bilan", type=Path)
    p.add_argument("destination", type=Path)
    a = p.parse_args()
    tracer(a.bilan, a.destination)

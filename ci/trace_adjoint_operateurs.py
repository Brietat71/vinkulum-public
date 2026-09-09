"""Temps et mémoire des produits adjoints, entre roues isolées."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter

LABELS={'ancien':('0.8.2 · matrices A/B','#454545'),
        'matrice':('0.9.0 · interface existante','#246ba7'),
        'produits':('0.9.0 · produits + état initial exact','#cc6825')}
TITLES={'pas':'Un produit de pas', 'partage':'Gradient · un EI partagé',
        'par_poutre':'Gradient · un EI par poutre'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bilan',type=Path);p.add_argument('sortie',type=Path)
    a=p.parse_args();report=json.loads(a.bilan.read_bytes())
    if not report['all_verified']:raise ValueError('mesures non vérifiées')
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for col,(family,title) in enumerate(TITLES.items()):
        for variant,(label,color) in LABELS.items():
            cs=sorted((c for c in report['configurations'] if (c['family'],c['variant'])==(family,variant)),key=lambda c:c['elements'])
            ns=[c['elements'] for c in cs];ts=[c['median_seconds'] for c in cs]
            axes[0,col].errorbar(ns,ts,yerr=[[t-min(c['seconds']) for t,c in zip(ts,cs)],
                [max(c['seconds'])-t for t,c in zip(ts,cs)]],label=label,color=color,marker='o',markersize=4,capsize=3)
            axes[1,col].plot(ns,[max(c['rss_kib'])/1024 for c in cs],color=color,marker='o',markersize=4)
        axes[0,col].set(title=title,ylabel='Temps médian du calcul (s)')
        axes[1,col].set(ylabel='Pic RSS maximal (Mio)',xlabel='Nombre de poutres')
        for ax in axes[:,col]:
            ax.set_xscale('log');ax.set_yscale('log');ax.grid(which='major',alpha=.2)
            ax.xaxis.set_major_locator(FixedLocator(ns));ax.xaxis.set_major_formatter(ScalarFormatter())
    axes[0,0].legend(fontsize=8,loc='upper left')
    fig.suptitle('Adjoint discret · 3 répétitions en processus frais · CPU 8, un fil\n'
                 'Chargement de SciPy inclus au premier recul ; gradients de 20 pas',fontsize=12)
    fig.savefig(str(a.sortie)+'.svg');fig.savefig(str(a.sortie)+'.png',dpi=150)
    svg=Path(str(a.sortie)+'.svg')
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')


if __name__=='__main__':main()

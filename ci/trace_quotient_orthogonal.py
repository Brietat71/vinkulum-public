"""Visualise l'archive du quotient : durées défavorables et remplissage inclus."""
import argparse
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('sortie', type=Path)
    args = parser.parse_args()
    data = json.loads((args.archive/'manifest.json').read_text())
    records = [r for r in data['essais'] if not r['echauffement']]
    cases = data['cas']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.hashsalt': 'vinkulum-quotient-2026'})
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout='constrained')
    labels, local, dense = [], [], []
    for family, size in cases:
        labels.append(f'{family.replace("_", " ")} · {size}')
        for variant, values in (('local', local), ('qr_dense', dense)):
            values.append(statistics.median(r['total_s'] for r in records if (r['famille'], r['taille'], r['variante']) == (family, size, variant)))
    y = list(range(len(cases)))
    axes[0].barh([i-.18 for i in y], local, height=.34, label='Quotient local Python', color='#236c83')
    axes[0].barh([i+.18 for i in y], dense, height=.34, label='QR dense des normales', color='#b35633')
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xscale('log')
    axes[0].set(title='Préparation + projection + réactions', xlabel='Secondes · médiane de 3 processus frais')
    axes[0].legend(loc='upper right', fontsize=9)
    for family, color, style in (('chaine', '#236c83', '-o'), ('cascade_tournee', '#b35633', '-s'), ('couplages', '#6c568c', '-^')):
        rows = [r for r in records if r['famille'] == family and r['variante'] == 'local' and r['repetition'] == 1]
        axes[1].loglog([r['forme'][1] for r in rows], [r['diagnostic']['nnz_facteur']/r['nnz'] for r in rows], style, color=color, label=family.replace('_', ' '))
    axes[1].set(title='Le facteur peut se remplir malgré un petit Q', xlabel='Nombre de coordonnées', ylabel='Coefficients du facteur R / coefficients de G')
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(True, axis='x' if ax is axes[0] else 'both', alpha=.2)
    fig.suptitle('Quotient orthogonal · sondage algébrique, aucun classement de solveurs', fontsize=13)
    args.sortie.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.sortie, metadata={'Date': None})
    if args.sortie.suffix == '.svg':
        args.sortie.write_text('\n'.join(line.rstrip() for line in args.sortie.read_text().splitlines())+'\n')
    fig.savefig(args.sortie.with_suffix('.png'), dpi=140)


if __name__ == '__main__':
    main()

"""Figure reproductible du coût à précision commune sur Princeton."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter

LABELS = {
    'vinkulum_milieu': ('Vinkulum · milieu', '#2875b9', '-'),
    'vinkulum_integree': ('Vinkulum · intégrée', '#178457', '-'),
    'mbdyn': ('MBDyn · beam3', '#292929', '-'),
    'exudyn_newton': ('Exudyn · Newton', '#d15b21', '-'),
    'exudyn_modifie': ('Exudyn · Newton modifié', '#a756a5', '-'),
    'exudyn_fast_newton': ('Exudyn fast · Newton', '#d15b21', '--'),
    'exudyn_fast_modifie': ('Exudyn fast · Newton modifié', '#a756a5', '--'),
}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bilan',type=Path); p.add_argument('sortie',type=Path)
    args=p.parse_args(); report=json.loads(args.bilan.read_text())
    fig,(ax,bars)=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.7,1]},layout='constrained')
    for variant,(label,color,style) in LABELS.items():
        cs=sorted((c for c in report['configurations'] if c['variant']==variant and 'median_seconds' in c),key=lambda c:c['n'])
        ax.loglog([1000*c['median_seconds'] for c in cs],[1e6*c['estimated_error_with_margin'] for c in cs],
                  style,marker='.',color=color,label=label,linewidth=1.4)
        best=report['best'][variant]
        if best:
            ax.plot(best['median_seconds']*1000,best['estimated_error_with_margin']*1e6,'*',color=color,markersize=12)
    ax.axhline(10,color='#666666',linestyle=':',linewidth=1)
    ax.set(xlabel='Temps total médian (ms)',ylabel='Erreur estimée au bout (µm)',title='Même seuil : 10 µm ; étoiles = réglages retenus')
    ax.grid(True,which='major',alpha=.18)
    ax.xaxis.set_major_locator(FixedLocator([20,50,100,200,500,1000]))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(fontsize=8,loc='upper right')
    families=['mbdyn','vinkulum','exudyn']
    for i,family in enumerate(families):
        candidates=[c for v,c in report['best'].items() if v.startswith(family) and c]
        if not candidates: continue
        c=min(candidates,key=lambda c:c['median_seconds']); ms=c['median_seconds']*1000
        label,color,_=LABELS[c['variant']]
        bars.barh(i,ms,color=color,height=.55)
        bars.errorbar(ms,i,xerr=[[ms-1000*min(c['wall_seconds'])],[1000*max(c['wall_seconds'])-ms]],color='black',capsize=3)
        bars.text(ms+6,i,f'{ms:.1f} ms',va='center',fontsize=10)
        bars.text(2,i+.36,label,fontsize=9)
    bars.set(yticks=[],xlabel='Temps total médian (ms)',title='Meilleur réglage testé par moteur',xlim=(0,340))
    bars.invert_yaxis(); bars.grid(axis='x',alpha=.18)
    fig.suptitle('Princeton statique · Vinkulum 0.8.2 / MBDyn / Exudyn 1.11.0',fontsize=13)
    fig.savefig(str(args.sortie)+'.svg')
    fig.savefig(str(args.sortie)+'.png',dpi=150)


if __name__=='__main__':
    main()

"""Courbes coût / erreur à partir du bilan mesuré, export SVG et PNG."""
import argparse
import json
from pathlib import Path


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bilan',type=Path); p.add_argument('sortie',type=Path)
    p.add_argument('--version',help='Version pour les archives anciennes dépourvues de cette métadonnée')
    args=p.parse_args(); report=json.loads(args.bilan.read_text())
    version=report['metadata'].get('vinkulum_version',args.version)
    if not version or (args.version and args.version!=version):
        p.error('version absente ou différente de celle des mesures')
    names=[n for n,c in report['cases'].items() if c['references_qualified']]
    if not names:
        p.error('aucune référence qualifiée à tracer')
    cols=min(3,len(names)); rows=(len(names)+cols-1)//cols
    fig,axes=plt.subplots(rows,cols,figsize=(5*cols,4.2*rows),squeeze=False)
    titles={'six_barres':'Mécanisme plan · 3 s','spatial':'Mécanisme spatial · 5 s','andrews':'Andrews · 20 ms', 'andrews_initialise':'Andrews · initialisation corrigée',
            'princeton':'Princeton · poutre milieu','princeton_integree':'Princeton · poutre intégrée'}
    for ax,name in zip(axes.flat,names):
        case=report['cases'][name]
        for engine,color in [('vinkulum','#1261a0'),('mbdyn','#d05a24')]:
            cs=[c for c in case['configurations'] if c['engine']==engine and c['candidate'] and c['successes']==c['attempts'] and c.get('estimated_error_with_margin',0)>0]
            cs=sorted({c['step']:c for c in cs}.values(),key=lambda c:c['median_seconds'])
            ax.loglog([c['median_seconds'] for c in cs],[c['estimated_error_with_margin'] for c in cs],'o-',label='Vinkulum' if engine=='vinkulum' else 'MBDyn',color=color)
            best=case[engine+'_best']
            if best:
                ax.scatter([best['median_seconds']],[best['estimated_error_with_margin']],marker='*',s=130,color=color,edgecolors='black',zorder=4)
        ax.axhline(case['spec']['target'],color='gray',ls='--',label='Seuil visé')
        ax.set(title=titles.get(name,name),xlabel='Temps total médian (s)',ylabel=f"Écart + marge des références ({case['spec']['unit']})")
        ax.xaxis.set_major_locator(LogLocator(base=10,subs=(1,2,5),numticks=12))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value,_:f'{value:g}'))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(True,which='both',alpha=.2); ax.legend(fontsize=8)
    for ax in list(axes.flat)[len(names):]:
        ax.set_visible(False)
    fig.suptitle(f'Vinkulum {version} / MBDyn installé — calculs séquentiels, un fil demandé\nÉtoiles : configurations retenues au seuil commun')
    fig.tight_layout()
    for extension in ('.svg','.png'): fig.savefig(Path(str(args.sortie)+extension),dpi=160)


if __name__=='__main__': main()

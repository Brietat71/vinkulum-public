"""Trace les dimensions et temps archivés, avec refus et échauffements explicites.

Les valeurs proviennent du manifeste et des fichiers bruts contrôlés par SHA256.
Le script ne lance aucune simulation. Les barres montrent les médianes ; les
moustaches montrent le minimum et le maximum des répétitions hors échauffement.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import textwrap


ROOT = Path(__file__).resolve().parents[1]
VARIANTES = ('krylov', 'hcb')
NOMS = {'krylov': 'Krylov adaptatif', 'hcb': 'HCB (réf. SciPy)'}
COULEURS = {'krylov': '#156f86', 'hcb': '#b45631'}


def lire_archive(archive):
    manifeste_path = archive/'manifest.json'
    if not manifeste_path.is_file():
        raise ValueError(f'Manifeste absent : {manifeste_path}')
    manifeste = json.loads(manifeste_path.read_text())
    cas = {c['nom']: c for c in manifeste['cas']}
    groupes = {(nom, variante): [] for nom in cas for variante in VARIANTES}
    vus = set()
    for essai in manifeste['essais']:
        cle = (essai['cas'], essai['variante'], essai['repetition'])
        if cle in vus or cle[:2] not in groupes:
            raise ValueError(f'Essai dupliqué ou inconnu : {cle}')
        vus.add(cle)
        fichier = archive/essai['fichier']
        contenu = fichier.read_bytes()
        if hashlib.sha256(contenu).hexdigest() != essai['sha256']:
            raise ValueError(f'Empreinte incorrecte : {fichier}')
        brut = json.loads(contenu)
        if brut['cas'] != cas[essai['cas']] or brut['variante'] != essai['variante']:
            raise ValueError(f'Identité incohérente : {fichier}')
        if essai['echauffement'] != (essai['repetition'] == 0):
            raise ValueError(f'Échauffement incohérent : {fichier}')
        for champ in ('total_s', 'preparation_s', 'requetes_s',
                      'reconstruction_s', 'borne_uniforme', 'erreur_schur_audit'):
            if not math.isfinite(brut[champ]) or brut[champ] < 0:
                raise ValueError(f'{champ} invalide : {fichier}')
        if not math.isclose(brut['total_s'], sum(brut[c] for c in
                            ('preparation_s', 'requetes_s', 'reconstruction_s')),
                            rel_tol=1e-12):
            raise ValueError(f'Durée totale incohérente : {fichier}')
        accepte = (brut['borne_uniforme'] <= brut['cas']['tolerance']
                   and brut['erreur_schur_audit'] <= brut['cas']['tolerance'])
        if brut['accepte'] != accepte:
            raise ValueError(f'Acceptation incohérente : {fichier}')
        if not essai['echauffement']:
            groupes[cle[:2]].append(brut)
    attendus = {(nom, variante, repetition) for nom in cas for variante in VARIANTES
                for repetition in range(manifeste['repetitions'])}
    if vus != attendus:
        raise ValueError('Archive incomplète : répétitions manquantes ou surnuméraires')
    nombres = {len(lignes) for lignes in groupes.values()}
    if len(nombres) != 1 or not next(iter(nombres)):
        raise ValueError('Répétitions hors échauffement incohérentes')
    return manifeste, groupes, next(iter(nombres))


def nombre(valeur, chiffres=3):
    return f'{valeur:.{chiffres}g}'.replace('.', ',')


def entier(valeur):
    return f'{valeur:,.0f}'.replace(',', '\u202f')


def famille(cas):
    return {'chaine': 'Chaîne', 'console': 'Console'}.get(cas['famille'], cas['famille'])


def etiquette(cas, lignes):
    bandes = {round(l['omega_max']/math.sqrt(l['lambda_min']), 12) for l in lignes}
    if len(bandes) != 1:
        raise ValueError(f'Bande incohérente : {cas["nom"]}')
    unite = 'masses' if cas['famille'] == 'chaine' else 'poutres'
    bande = nombre(next(iter(bandes)), 3)
    return f'{famille(cas)}\n{entier(cas["n"])} {unite}\nΩ/√λ* = {bande}'


def statistiques(lignes, champ, facteur=1.):
    valeurs = [l[champ]*facteur for l in lignes]
    mediane = statistics.median(valeurs)
    return mediane, mediane-min(valeurs), max(valeurs)-mediane


def texte_refus(cas, groupes):
    comptes = []
    for variante in VARIANTES:
        lignes = groupes[(cas['nom'], variante)]
        refuses = sum(not l['accepte'] for l in lignes)
        comptes.append(f'{NOMS[variante]} : {refuses}/{len(lignes)} refusés')
    return (f'{famille(cas)} {entier(cas["n"])} à tolérance {nombre(cas["tolerance"])} '
            f'— {" ; ".join(comptes)}. Cas hors classement.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', nargs='?', type=Path,
                        default=ROOT/'docs/bancs/ports-krylov-2026')
    parser.add_argument('sortie', nargs='?', type=Path,
                        default=ROOT/'docs/figures/ports-krylov-2026',
                        help='préfixe des sorties .svg et .png')
    args = parser.parse_args()
    try:
        manifeste, groupes, repetitions = lire_archive(args.archive)
        classes, ecartes = [], []
        for cas in manifeste['cas']:
            if cas['nom'].endswith('_arrondi'):
                ecartes.append(cas)
                continue
            lignes = [l for v in VARIANTES for l in groupes[(cas['nom'], v)]]
            if not all(l['accepte'] for l in lignes):
                raise ValueError(f'Cas non accepté hors série arrondi : {cas["nom"]}')
            if any(l['directions'] <= 0 or l['total_s'] <= 0 for l in lignes):
                raise ValueError(f'Valeur non positive sur une échelle log : {cas["nom"]}')
            classes.append(cas)
        if not classes:
            raise ValueError('Aucun cas accepté à tracer')
        mesures = [l for c in classes for v in VARIANTES for l in groupes[(c['nom'], v)]]
        tolerances = {l['cas']['tolerance'] for l in mesures}
        requetes = {l['nombre_requetes'] for l in mesures}
        if len(tolerances) != 1 or len(requetes) != 1:
            raise ValueError('Tolérance ou nombre de requêtes différents entre cas classés')
        tolerance, requetes = next(iter(tolerances)), next(iter(requetes))
        labels = [etiquette(c, [l for v in VARIANTES for l in groupes[(c['nom'], v)]])
                  for c in classes]
    except (ValueError, KeyError, OSError) as exc:
        parser.error(str(exc))

    # Import après validation : une archive absente ne déclenche aucun tracé.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, LogLocator

    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 11,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.titleweight': 'bold', 'axes.labelcolor': '#27343b',
        'text.color': '#27343b', 'axes.edgecolor': '#bdc6ca',
        'svg.hashsalt': 'vinkulum-ports-krylov-2026', 'svg.fonttype': 'none',
    })
    fig, axes = plt.subplots(1, 2, figsize=(16, 8.5))
    fig.subplots_adjust(left=.065, right=.985, bottom=.31, top=.76, wspace=.19)
    fig.suptitle('Réduction par ports à tolérance commune',
                 x=.065, y=.965, ha='left', fontsize=20, fontweight='bold')
    fig.text(.065, .913,
             f'Tolérance du Schur normalisé : {nombre(tolerance)} · '
             f'{entier(requetes)} requêtes dans la bande + construction + borne + une reconstruction de champ',
             ha='left', fontsize=11)
    fig.text(.065, .877,
             f'Médianes de {repetitions} répétitions hors échauffement ; moustaches = minimum–maximum. '
             'Bande relative indiquée sous chaque cas.', ha='left', fontsize=10.5)

    x = list(range(len(classes)))
    for ax, champ, facteur, titre, unite in (
            (axes[0], 'directions', 1., 'Dimension intérieure retenue', 'Directions intérieures'),
            (axes[1], 'total_s', 1000., 'Coût total par problème', 'Temps total (ms)')):
        for variante, decalage in zip(VARIANTES, (-.18, .18)):
            stats = [statistiques(groupes[(c['nom'], variante)], champ, facteur) for c in classes]
            medianes = [s[0] for s in stats]
            bas = [s[1] for s in stats]
            haut = [s[2] for s in stats]
            positions = [xi+decalage for xi in x]
            ax.bar(positions, medianes, width=.32, color=COULEURS[variante],
                   label=NOMS[variante], zorder=3)
            ax.errorbar(positions, medianes, yerr=[bas, haut], fmt='none',
                        ecolor='#27343b', capsize=3, linewidth=1., zorder=4)
            for pos, mediane, marge in zip(positions, medianes, haut):
                texte = entier(mediane) if champ == 'directions' else nombre(mediane)
                ax.annotate(texte, (pos, mediane+marge), xytext=(0, 6),
                            textcoords='offset points', ha='center', va='bottom', fontsize=9)
        ax.set_yscale('log')
        ax.set_xticks(x, labels, fontsize=9.5)
        ax.tick_params(axis='x', length=0, pad=10)
        ax.set_title(titre, loc='left', fontsize=14, pad=15)
        ax.set_ylabel(unite)
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda valeur, _: nombre(valeur)))
        ax.grid(True, axis='y', which='major', color='#d9dfe2', linewidth=.8, zorder=0)
        ax.set_axisbelow(True)
        ax.margins(y=.26, x=.04)
    handles, noms = axes[0].get_legend_handles_labels()
    fig.legend(handles, noms, loc='upper right', bbox_to_anchor=(.987, .852),
               ncol=len(VARIANTES), frameon=False, fontsize=11)

    notes_refus = ' '.join(texte_refus(c, groupes) for c in ecartes)
    if notes_refus:
        fig.text(.065, .177, '\n'.join(textwrap.wrap(notes_refus, width=170)),
                 fontsize=10.5, va='top', color='#843f28')
    fig.text(.065, .105,
             'Imports, assemblage et audit indépendant exclus des durées. '
             'Nombre de modes HCB présélectionné, coût de recherche exclu.', fontsize=10)
    fig.text(.065, .073,
             'Référence HCB exécutée avec SciPy ; aucun temps attribué à Exudyn, MBDyn ou Simpack. '
             'Bornes évaluées sans certification des arrondis.', fontsize=10)
    try:
        source = args.archive.resolve().relative_to(ROOT)
    except ValueError:
        source = args.archive
    fig.text(.065, .037, f'Source : {source}/manifest.json et essais bruts associés.',
             fontsize=9, color='#64757e')
    prefixe = args.sortie.with_suffix('') if args.sortie.suffix in ('.svg', '.png') else args.sortie
    prefixe.parent.mkdir(parents=True, exist_ok=True)
    svg, png = prefixe.with_suffix('.svg'), prefixe.with_suffix('.png')
    empreinte = hashlib.sha256((args.archive/'manifest.json').read_bytes()).hexdigest()
    description = f'Archive {source}, manifeste SHA256 {empreinte}. Médianes hors échauffement.'
    fig.savefig(svg, metadata={'Date': None, 'Title': 'Vinkulum : réduction par ports',
                              'Description': description, 'Creator': 'Vinkulum'})
    svg.write_text('\n'.join(ligne.rstrip() for ligne in svg.read_text().splitlines())+'\n')
    fig.savefig(png, dpi=180, metadata={'Description': description, 'Software': 'Vinkulum'})
    plt.close(fig)
    print(svg)
    print(png)


if __name__ == '__main__':
    main()

"""Requalification des assemblages et rejeu des neuf paires de certificats."""
import io,json,tempfile,math
from pathlib import Path
import numpy as np
from archive_masse_comparee import archiver,digest
from archive_inertie_binaire import inventaire,octets,lire_json,preuve_stable
from experience_assemblage_complements import SOURCES,CI,VARIANTES,PARAMETRES
from experience_controle_facteurs import qualifier,qualifier_controle,THREADS,FILS_SUPPLEMENTAIRES,BUDGETS_INERTIE
from modeles_assemblage_complements import charger_cas
from comparaison_masse import compiler_comparaison,BibliothequeComparaison,certifier_comparaison_masse
from inertie_binaire_compile import BibliothequeInertie,certifier_inertie_compilee

ARCHIVE=CI.parent/'docs/bancs/assemblage-complements-2026'

def verifier_rapport(rows,passages=4):
    plan=[(p,n,v) for p in range(passages) for n in (32,128,512) for v in VARIANTES]
    if [(r['passage'],r['n'],r['variante']) for r in rows]!=plan:raise ValueError('plan incomplet ou dupliqué')
    envs=set();admitted=controls=0;sources=rows[0]['sources_sha256'];binary=rows[0]['bibliotheque']['sha256']
    for r in rows:
        if r['role']!=('chauffe' if r['passage']==0 else 'mesure') or r['parametres']!=PARAMETRES:raise ValueError('rôle ou paramètres différents')
        e=r['environnement']
        if e['cpu']!=[8] or e['fils']!=THREADS or e['fils_supplementaires']!=FILS_SUPPLEMENTAIRES:raise ValueError('budget de calcul différent')
        envs.add(json.dumps({k:e[k] for k in ('python','numpy','scipy','executable','distributions','bibliotheques_numeriques_sha256')},sort_keys=True))
        if r['sources_sha256']!=sources or r['bibliotheque']['sha256']!=binary:raise ValueError('sources ou binaire instables')
        if r['certification_machine_reponses'] is not False:raise ValueError('certification non justifiée')
        q=qualifier(r)
        if q!=r['qualifie_champs']:raise ValueError('qualification du champ incohérente')
        admitted+=q
        if r['variante']=='assemblage_controle':
            qc=qualifier_controle(r)
            if list(qc)!=r['qualifie_controle']:raise ValueError('qualification du contrôle incohérente')
            controls+=qc[0]
            if r['statut']=='termine' and (len(r['certificats'])!=3 or r['taille_conservee']!=9):raise ValueError('couverture de pièces incohérente')
    if len(envs)!=1:raise ValueError('environnements numériques différents')
    return dict(essais=len(rows),champs_admis=admitted,controles_admis=controls)


def verifier_archive(racine=ARCHIVE,*,rejouer=True):
    racine=Path(racine);manifest=inventaire(racine)
    for group in ('campagne','pilote','references','conformite'):
        for name,expected in lire_json(racine,group+'/manifest.json').items():
            full=group+'/'+name
            actual=manifest['exclus_npy'][full]['source_sha256'] if full in manifest['exclus_npy'] else digest(octets(racine,full))
            if actual!=expected:raise ValueError('objet différent du manifeste original')
    rows=lire_json(racine,'campagne/bilan.json');result=verifier_rapport(rows)
    audit=lire_json(racine,'conformite/audit.json')
    if audit['certification_machine'] is not False or len(audit['essais'])!=len(rows):raise ValueError('audit de conformité incomplet')
    for a,r in zip(audit['essais'],rows,strict=True):
        if a['dossier']!=r['dossier'] or a['champs_sha256']!=r['champs_sha256']:raise ValueError('conformité jugée sur un autre champ')
        if len(a['colonnes'])!=257 or len(a['operateurs'])!=257 or any(len(row)!=6 for row in a['colonnes']):raise ValueError('couverture de conformité incomplète')
        values=[v for row in a['colonnes'] for v in row]
        if not all(math.isfinite(v) and v>=0 for v in values+a['operateurs']):raise ValueError('conformité non finie')
        if max(values)!=a['maximum_colonnes'] or max(a['operateurs'])!=a['maximum_operateur']:raise ValueError('maximum de conformité incohérent')
    name='audit_jonctions_complements.py'
    if digest(octets(racine,'conformite/'+name))!=digest((CI/name).read_bytes()):raise ValueError('source de conformité différente')
    verifier_rapport(lire_json(racine,'pilote/bilan.json'),passages=1)
    identity=lire_json(racine,'campagne/build/identite.json')
    if digest(octets(racine,'campagne/build/comparaison_masse.so'))!=identity['bibliotheque_sha256']:raise ValueError('binaire capturé altéré')
    for r in rows:
        prefix='campagne/'+r['dossier'];single=lire_json(racine,prefix+'/resultat.json')
        if single!={k:v for k,v in r.items() if k not in ('passage','role','dossier')}:raise ValueError('résultat individuel différent')
        if lire_json(racine,prefix+'.terminal.json')['code']!=0:raise ValueError('processus en échec')
        if manifest['exclus_npy'][prefix+'/champs.npy']['source_sha256']!=r['champs_sha256']:raise ValueError('champ exclu incohérent')
        for name,expected in r['sources_sha256'].items():
            if digest(octets(racine,'campagne/sources/'+name))!=expected or digest((CI/name).read_bytes())!=expected:raise ValueError('source différente de la mesure')
        case=f'references/n{r["n"]}'
        if lire_json(racine,case+'/modele.json')!=r['modele']:raise ValueError('modèle différent')
        for name,expected in r['modele']['fichiers_sha256'].items():
            if digest(octets(racine,case+'/'+name))!=expected:raise ValueError('branche ou jonction différente')
        ref=lire_json(racine,case+'/reference.json')
        if (ref['precisions']!=[70,90] or ref['certification_machine'] is not False
            or not 0<=ref['ecart_relatif_max']<=1e-8
            or ref['ecart_relatif_max']!=max(*ref['maxima'].values(),*ref['maxima_operateurs'].values())):raise ValueError('référence non qualifiée')
        if r['reference_sha256']!=ref['reference_sha256'] or manifest['exclus_npy'][case+'/reference90.npy']['source_sha256']!=r['reference_sha256']:raise ValueError('référence différente')
        if r['bibliotheque']['sha256']!=identity['bibliotheque_sha256']:raise ValueError('bibliothèque différente')
    cache={}
    if rejouer:
        with tempfile.TemporaryDirectory(prefix='vinkulum-rejeu-assemblage-') as tmp:
            lib=compiler_comparaison(Path(tmp)/'build');bib=BibliothequeComparaison(lib);inertia=BibliothequeInertie(lib)
            for r in rows:
                if r['variante']!='assemblage_controle' or r['statut']!='termine':continue
                n=r['n'];pieces,j=charger_cas(racine/f'references/n{n}')
                for index,certs in enumerate(r['certificats']):
                    raw=octets(racine,'campagne/'+r['dossier']+f'/base{index}.npz')
                    with np.load(io.BytesIO(raw),allow_pickle=False) as z:b=z['B']
                    key=(n,index,digest(b.tobytes()))
                    if key not in cache:
                        d,m,h=pieces[index];mi=m[:-6,:-6].tocsr()
                        comp=certifier_comparaison_masse(mi,.49,1.51,bibliotheque=bib)
                        actual=certifier_inertie_compilee(d[:,:-6],mi,b,float((2*np.pi*80)**2),bibliotheque=inertia,
                            permutation=certs['inertie']['permutation_physique'],**BUDGETS_INERTIE)
                        cache[key]=(preuve_stable(comp),preuve_stable(actual))
                    if cache[key]!=(preuve_stable(certs['comparaison']),preuve_stable(certs['inertie'])):raise ValueError('preuve différente du rejeu')
    return dict(**result,paires_certificats_rejouees=len(cache))

if __name__=='__main__':print(verifier_archive())

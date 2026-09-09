"""Archive et rejoue l'attribution discrète, sans la confondre avec une preuve ODE."""
import argparse
import gzip
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

H=('0.00100','0.00050','0.00025')
FRAMES=('spatial','materiel','combine')
CASES=tuple(f'rotation-{frame}-{h}' for h in H for frame in FRAMES)
TRACES=tuple(f'{frame}-{h}.json.gz' for h in H for frame in ('initial',*FRAMES))


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def rational(x):
    if not isinstance(x,dict) or set(x)!={'n','d'}:raise ValueError('rationnel absent')
    n,d=int(x['n']),int(x['d'])
    if d<=0:raise ValueError('dénominateur non positif')
    return F(n,d)


def require(ok,message):
    if not ok:raise ValueError(message)


def bilan(directory):
    """Résumé dérivé des rationnels archivés, sans nouveau critère de succès."""
    directory=Path(directory);rows=[]
    for name in CASES:
        result=read(directory/'resultats'/(name+'.json'));rv=result['rotation_vitesse']
        omega=rv['pics_tous_pas']['omega'];contributions=omega['contributions']
        lower=rational(contributions[2][0])
        others=sum(rational(contributions[i][1]) for i in (0,1,3))
        rows.append(dict(cadre=result['cadre'],pas=result['pas_demande'],pas_verifies=result['pas_verifies'],
                         rotation_max=float(rational(rv['pics_tous_pas']['rotation']['norme'][1])),
                         omega_max=float(rational(omega['norme'][1])),
                         contributions_au_pic_omega=[float(rational(v[1])) for v in contributions],
                         defauts_dominants_au_pic_omega=lower>others,
                         largeur_moment=float(rational(result['largeur_max'])),
                         largeur_rotation=float(rational(rv['largeur_rotation_max'])),
                         largeur_omega=float(rational(rv['largeur_omega_max']))))
    return dict(schema='vinkulum.bilan.attribution_reperes_em.1',cas=rows,
                valeurs_flottantes='approximations affichées ; les bornes normatives sont les rationnels des résultats',
                pas_comparaisons=sum(row['pas_verifies'] for row in rows),
                defauts_dominants_aux_neuf_pics_omega=all(row['defauts_dominants_au_pic_omega'] for row in rows),
                definition_dominance='borne inférieure de la norme de la contribution des défauts supérieure à la somme des bornes supérieures des trois autres contributions, au pic de vitesse de chaque comparaison',
                borne_erreur_ode_continue=False)


def judge(result,name,trace_hashes,source_hashes):
    _,frame,h=name.split('-')
    require(result['cadre']==frame and result['pas_demande']==float(h),'cas mal identifié')
    expected=round(20./float(h))
    require(result['complet'] is True and result['pas_total']==result['pas_verifies']==expected,'horizon incomplet')
    require(result['sources_sha256']==[trace_hashes[f'initial-{h}.json.gz'],trace_hashes[f'{frame}-{h}.json.gz']],'traces non rattachées')
    require(result['programme_sha256']==source_hashes['ci/attribue_moments_em.py'],'programme moment non rattaché')
    budget=F(64,2**52)
    require(rational(result['budget_local'])==budget,'budget local changé')
    require(rational(result['budget_largeur'])==F(1,10**20),'budget moment changé')
    require(result['precision_bits']==160,'précision changée')
    require(all(0<=rational(x)<=budget for x in result['defauts_locaux_max']),'défauts de moment excessifs')
    require(len(result['defauts_locaux_max'])==2,'défauts absents')
    require(0<=rational(result['largeur_max'])<=F(1,10**20),'enclosure moment excessive')
    rv=result['rotation_vitesse']
    require(rv['programme_sha256']==source_hashes['ci/attribue_rotation_em.py'],'programme rotation non rattaché')
    require(rv['pas_verifies']==expected,'rotation incomplète')
    require(rv['borne_erreur_ode_continue'] is False,'portée ODE abusive')
    require(rv['determinants_positifs'] is True,'orientation refusée')
    require(rational(rv['budget_local'])==budget,'budget rotation local changé')
    for field in ('defaut_rotation_local_max','defaut_omega_local_max','orthogonalite_native_max','orthogonalite_reperes_max'):
        require(0<=rational(rv[field])<=budget,'critère local refusé: '+field)
    for field,den in (('rotation',20),('omega',18)):
        require(rational(rv['budget_'+field])==F(1,10**den),'budget enclosure changé')
        require(0<=rational(rv['largeur_'+field+'_max'])<=F(1,10**den),'enclosure excessive')
    require(result['attribution_moment_verifiee'] is True and result['attribution_rotation_vitesse_verifiee'] is True and rv['verifie'] is True,'attribution refusée')


def construire(source,destination):
    source,destination=Path(source),Path(destination)
    require(not destination.exists(),'destination déjà existante')
    provenance=read(source/'provenance-campagne.json')
    hashes={name:sha(source/'traces'/name) for name in TRACES}
    for name in CASES:judge(read(source/(name+'.json')),name,hashes,provenance['sources'])
    for name,digest in provenance['sources'].items():
        require(sha(source/'sources-campagne'/name)==digest,'snapshot altéré: '+name)
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
        out=Path(temporary)/'archive';out.mkdir()
        shutil.copytree(source/'sources-campagne',out/'sources')
        (out/'traces').mkdir();(out/'resultats').mkdir()
        for name in TRACES:shutil.copyfile(source/'traces'/name,out/'traces'/name)
        for name in CASES:
            for suffix in ('.json','.log'):shutil.copyfile(source/(name+suffix),out/'resultats'/(name+suffix))
        for original,target in [('provenance-campagne.json','provenance.json'),('traces/manifest.json','capture.json'),('baseline/manifest.json','baseline.json'),('rattachement-historique.json','historique.json'),('candidate-build.log','construction-candidate.log'),('tests-package.log','tests-package.log')]:
            shutil.copyfile(source/original,out/target)
        (out/'bilan.json').write_text(json.dumps(bilan(out),indent=2)+'\n')
        files={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
        manifest=dict(schema='vinkulum.archive.attribution_reperes_em.1',version_candidate='0.18.0',
                      cas=list(CASES),traces=list(TRACES),fichiers_sha256=files,
                      attribution_discrete_verifiee=True,borne_erreur_ode_continue=False)
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        verifier(out)
        out.rename(destination)
    return manifest


def verifier(directory):
    directory=Path(directory);manifest=read(directory/'manifest.json')
    require(manifest['schema']=='vinkulum.archive.attribution_reperes_em.1','schéma inconnu')
    require(manifest['cas']==list(CASES) and manifest['traces']==list(TRACES),'inventaire incomplet')
    expected={str(p.relative_to(directory)) for p in directory.rglob('*') if p.is_file() and p!=directory/'manifest.json'}
    require(set(manifest['fichiers_sha256'])==expected,'inventaire de fichiers différent')
    for name,digest in manifest['fichiers_sha256'].items():
        require(sha(directory/name)==digest,'fichier altéré: '+name)
    provenance=read(directory/'provenance.json')
    for name,digest in provenance['sources'].items():require(sha(directory/'sources'/name)==digest,'source altérée')
    hashes={name:sha(directory/'traces'/name) for name in TRACES}
    capture=read(directory/'capture.json');baseline=read(directory/'baseline.json')
    require(capture['version']=='0.18.0' and baseline['version']=='0.17.0','roues mal identifiées')
    require(capture['extension_sha256']==provenance['candidate_extension_sha256'],'extension non rattachée')
    require(capture['programme_sha256']==provenance['sources']['ci/trace_moments_em.py'],'capture non rattachée au programme')
    require(capture['modele_sha256']==provenance['sources']['ci/diagnostic_rotation.py'],'modèle non rattaché')
    require({r['fichier'] for r in capture['cas']}==set(TRACES) and len(capture['cas'])==12,'capture incomplète')
    require({r['fichier'] for r in baseline['cas']}==set(TRACES) and len(baseline['cas'])==12,'baseline incomplète')
    for row in capture['cas']:
        require(row['sha256']==hashes[row['fichier']],'capture non rattachée')
        for key in ('comparaison_meme_roue','comparaison_roue_publiee'):
            require(row[key]['composantes_differentes_bits']==0 and row[key]['ecart_absolu_max']==0,'instrumentation différente')
    historical=read(directory/'historique.json')
    require(historical['toutes_identiques'] is True and len(historical['comparaisons'])==12,'rattachement historique incomplet')
    require(all(r['identiques'] is True and r['historique']==r['trace'] and r['differences']==[0.,0.] for r in historical['comparaisons']),'observations historiques différentes')
    for name in CASES:judge(read(directory/'resultats'/(name+'.json')),name,hashes,provenance['sources'])
    require(read(directory/'bilan.json')==bilan(directory),'bilan différent des résultats')
    require(manifest['attribution_discrete_verifiee'] is True and manifest['borne_erreur_ode_continue'] is False,'portée incorrecte')
    return dict(cas=len(CASES),traces=len(TRACES),controle='intégrité et décisions archivées ; pas un rejeu rationnel')


def verifier_observations(directory):
    """Relit les points aux pics ; ne remplace pas les identités de chaque pas."""
    from observations_reperes_archivees import ecarts
    from attribue_rotation_em import norm_bounds,sqrt_bounds
    from attribue_moments_em import iadd
    directory=Path(directory);verifier(directory)
    count=0
    for h in H:
        def load(frame):
            return json.loads(gzip.decompress((directory/'traces'/f'{frame}-{h}.json.gz').read_bytes()))
        base=load('initial')['trace']['echantillons']
        for frame in FRAMES:
            data=load(frame);rows=data['trace']['echantillons']
            rv=read(directory/'resultats'/f'rotation-{frame}-{h}.json')['rotation_vitesse']
            for mode in ('pics_tous_pas','pics_echantillonnage_historique'):
                for field in ('rotation','omega'):
                    peak=rv[mode][field];index=peak['pas']
                    require(type(index)==int and 0<=index<len(rows),'indice de pic invalide')
                    a,b=base[index],rows[index]
                    dr,dw=ecarts(a,b,data['monde'],data['matiere'])
                    observed=dr if field=='rotation' else dw
                    exact=[F(v) for v in observed]
                    square=sum(v*v for v in exact)
                    require(rational(peak['carre'])==square,'pic différent de la trace')
                    require(tuple(map(rational,peak['norme']))==sqrt_bounds(square),'norme de pic incorrecte')
                    vectors=[[(rational(lo),rational(hi)) for lo,hi in v] for v in peak['vecteurs']]
                    require(len(vectors)==4 and all(len(v)==len(exact) for v in vectors),'contributions incomplètes')
                    require(all(lo<=hi for v in vectors for lo,hi in v),'intervalle inversé')
                    total=vectors[0]
                    for vector in vectors[1:]:total=iadd(total,vector)
                    require(all(lo<=v<=hi for (lo,hi),v in zip(total,exact)),'observation hors contributions')
                    require(len(peak['contributions'])==4,'normes de contribution absentes')
                    for vector,bounds in zip(vectors,peak['contributions']):
                        require(tuple(map(rational,bounds))==norm_bounds(vector),'norme de contribution incorrecte')
                    count+=1
    return dict(pics_rattaches=count,rejeu_complet=False)


def rejouer_cas(directory,name,out,limite=None):
    directory=Path(directory).resolve()
    require(name in CASES,'cas inconnu')
    _,frame,h=name.split('-')
    command=[sys.executable,'-B',str(directory/'sources/ci/attribue_rotation_em.py'),
             str(directory/'traces'/f'initial-{h}.json.gz'),str(directory/'traces'/f'{frame}-{h}.json.gz'),
             '--sortie',str(out)]
    if limite is not None:command+=['--limite',str(limite)]
    # Les imports du rejeu ne doivent pas ajouter de __pycache__ à l'archive.
    subprocess.run(command,check=True)
    return read(out)


def rejouer(directory):
    directory=Path(directory).resolve();verifier(directory)
    with tempfile.TemporaryDirectory() as temporary:
        for name in CASES:
            out=Path(temporary)/(name+'.json')
            new=rejouer_cas(directory,name,out);old=read(directory/'resultats'/(name+'.json'))
            require(new==old,'rejeu différent: '+name)
    return dict(cas_rejoues=len(CASES),identiques=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--source');group.add_argument('--verifier');group.add_argument('--rejouer')
    p.add_argument('--sortie');p.add_argument('--observations',action='store_true');a=p.parse_args()
    if a.source:
        if not a.sortie:p.error('--sortie requis avec --source')
        construire(a.source,a.sortie);print(verifier(a.sortie))
    elif a.verifier:
        print(verifier(a.verifier))
        if a.observations:print(verifier_observations(a.verifier))
    else:print(rejouer(a.rejouer))

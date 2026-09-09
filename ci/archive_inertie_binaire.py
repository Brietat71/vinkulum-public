"""Archives compactes : requalification et rejeu des six certificats connus.

Aucune source ni bibliothèque capturée n'est exécutée. Le C++ courant du
projet est compilé dans un dossier temporaire. Les grands champs exclus
ne sont pas rejugés ; leurs identités et diagnostics restent vérifiables.
"""
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile

import numpy as np

from experience_inertie_binaire import (analyser,plan,sources,qualifier,qualifier_controle,
                                       GAMMA_INERTIE,BUDGETS_INERTIE)
from inertie_binaire_compile import compiler,BibliothequeInertie,certifier_inertie_compilee
from inertie_complement_dirigee import certifier_inertie_complement
from ordre_separateurs import ordre_separateurs
from confronte_ports_exudyn import lire

ARCHIVE=Path(__file__).resolve().parent.parent/'docs/bancs/inertie-binaire-2026'
SONDES=ARCHIVE.with_name('inertie-binaire-sondes-2026')


def sha(data):return hashlib.sha256(data).hexdigest()


def octets(racine,nom):
    p=Path(nom)
    if p.is_absolute() or '..' in p.parts:raise ValueError('chemin relatif requis')
    compressed=racine/(nom+'.gz')
    return gzip.decompress(compressed.read_bytes()) if compressed.is_file() else (racine/nom).read_bytes()


def lire_json(racine,nom):return json.loads(octets(racine,nom))


def inventaire(racine):
    racine=Path(racine)
    manifest=json.loads((racine/'manifest.json').read_text())
    if manifest['schema']!=1:raise ValueError('schéma de manifeste inattendu')
    present={str(p.relative_to(racine)) for p in racine.rglob('*') if p.is_file() and p!=racine/'manifest.json'}
    if present!=set(manifest['fichiers']):raise ValueError('inventaire incomplet ou supplémentaire')
    for name,metadata in manifest['fichiers'].items():
        if Path(name).is_absolute() or '..' in Path(name).parts:raise ValueError('chemin relatif requis')
        data=(racine/name).read_bytes()
        if sha(data)!=metadata['sha256'] or len(data)!=metadata['octets']:
            raise ValueError('octets altérés : '+name)
        if 'source_sha256' in metadata:
            if metadata['compression'] not in (None,'gzip'):raise ValueError('compression inconnue')
            raw=gzip.decompress(data) if metadata['compression']=='gzip' else data
            if sha(raw)!=metadata['source_sha256'] or len(raw)!=metadata['source_octets']:
                raise ValueError('copie différente de ses octets d’origine')
    return manifest


def preuve_stable(cert):
    return {k:v for k,v in cert.items() if k not in ('preparation_s','phases_s','phases_natives_s','bibliotheque')}


def verifier_certificats(rapport,racine,bibliotheque):
    cache={}
    for essai in rapport['essais']:
        result=essai['resultat']
        if result['variante']=='lu_corrigee' or result['statut']!='termine':continue
        n=result['n'];cert=result['certificat_complement'];meta=result['bases']
        raw=octets(racine,'essais/'+essai['nom']+'/'+meta['fichier'])
        if sha(raw)!=meta['sha256']:raise ValueError('base physique altérée')
        with np.load(io.BytesIO(raw),allow_pickle=False) as z:b,phi=z['B'],z['Phi']
        if b.shape!=(6*n-6,1) or phi.shape!=b.shape or not np.all(np.isfinite(b)) or not np.all(np.isfinite(phi)):
            raise ValueError('forme ou finitude B/Phi invalide')
        bsha=sha(np.asarray(b,dtype='<f8').tobytes());psha=sha(np.asarray(phi,dtype='<f8').tobytes())
        if bsha!=meta['b_sha256'] or bsha!=cert['contraintes_sha256'] or psha!=meta['phi_sha256']:
            raise ValueError('B ou Phi différents du certificat')
        key=(n,result['variante'],bsha)
        if key not in cache:
            d,m,*_=lire(racine/'entrees'/f'n{n}-f40.npz');di=d[:,:6*n-6];mi=m[:6*n-6,:6*n-6].tocsr()
            if result['variante']=='binaire_controle':
                perm=ordre_separateurs(di,mi,taille_feuille=8)
                actual=certifier_inertie_compilee(di,mi,b,GAMMA_INERTIE,bibliotheque=bibliotheque,
                                                permutation=perm,**BUDGETS_INERTIE)
            else:actual=certifier_inertie_complement(di,mi,b,GAMMA_INERTIE,precision=32,**BUDGETS_INERTIE)
            cache[key]=preuve_stable(actual)
        if preuve_stable(cert)!=cache[key]:raise ValueError('preuve différente du rejeu local')
    return len(cache)


def verifier_rapport(rapport):
    bilan=analyser(rapport)
    if rapport['statut']!='termine' or bilan!=rapport['analyse']:raise ValueError('bilan incohérent')
    if [(x['resultat']['n'],x['resultat']['variante'],x['passage']) for x in rapport['essais']]!=plan():
        raise ValueError('ordre des essais incohérent')
    for n in (32,128,512):
        rows=[x['resultat'] for x in rapport['essais'] if x['resultat']['n']==n and x['resultat']['variante']!='lu_corrigee']
        if len({x['champs']['sha256'] for x in rows})!=1:raise ValueError('champs différents entre contrôleurs')
        for key in ('normes_candidates','retours_par_frequence'):
            if any(x[key]!=rows[0][key] for x in rows):raise ValueError('normes ou majorants différents entre contrôleurs')
    return bilan


def verifier_archive(racine=ARCHIVE,*,rejouer=True):
    racine=Path(racine);manifest=inventaire(racine);rapport=lire_json(racine,'rapport.json')
    original=lire_json(racine,'manifest.json')['fichiers_sha256']
    for name,expected in original.items():
        if name in manifest['exclus_npy']:
            if manifest['exclus_npy'][name]['source_sha256']!=expected:raise ValueError('NPY exclu différent du manifeste original')
        elif sha(octets(racine,name))!=expected:raise ValueError('objet différent du manifeste original')
    if rapport['sources_sha256']!=sources():raise ValueError('sources courantes différentes des sources mesurées')
    for name,expected in rapport['sources_sha256'].items():
        if sha(octets(racine,'sources/'+name))!=expected:raise ValueError('source capturée altérée')
    comp=lire_json(racine,'build/identite.json')
    if (comp!=rapport['compilation'] or sha(octets(racine,'build/inertie_binaire.so'))!=comp['bibliotheque_sha256']
            or sha(octets(racine,'build/inertie_binaire_native.cpp'))!=comp['source_sha256']
            or comp['source_sha256']!=rapport['sources_sha256']['inertie_binaire_native.cpp']):
        raise ValueError('compilation incohérente')
    bilan=verifier_rapport(rapport)
    if bilan!=lire_json(racine,'bilan.json'):raise ValueError('bilan séparé incohérent')
    exclus_attendus=set()
    for e in rapport['essais']:
        prefix='essais/'+e['nom']+'/';r=e['resultat'];raw=octets(racine,prefix+'resultat.json')
        if sha(raw)!=e['resultat_sha256'] or json.loads(raw)!=r:raise ValueError('résultat individuel incohérent')
        if sha(octets(racine,'essais/'+e['nom']+'.log'))!=e['journal_sha256']:raise ValueError('journal altéré')
        terminal=lire_json(racine,'essais/'+e['nom']+'.terminal.json')
        if terminal!=e['terminal'] or terminal['retour']!=0 or terminal['attente_terminee'] is not True:
            raise ValueError('processus non terminé correctement')
        name=prefix+r['champs']['fichier'];exclus_attendus.add(name)
        if manifest['exclus_npy'][name]['source_sha256']!=r['champs']['sha256']:raise ValueError('champ exclu incohérent')
        if 'resultat_historique' in r:
            h=r['resultat_historique']
            if sha(octets(racine,prefix+h['fichier']))!=h['sha256']:raise ValueError('témoin brut altéré')
    if set(manifest['exclus_npy'])!=exclus_attendus:raise ValueError('exclusions NPY inattendues')
    replays=0
    if rejouer:
        with tempfile.TemporaryDirectory(prefix='vinkulum-rejeu-inertie-binaire-') as tmp:
            bib=BibliothequeInertie(compiler(Path(tmp)/'build'))
            replays=verifier_certificats(rapport,racine,bib)
    return dict(essais=len(rapport['essais']),admis=sum(qualifier(x['resultat']) for x in rapport['essais']),
        controles_admis=sum(qualifier_controle(x['resultat'])[0] for x in rapport['essais']
                            if x['resultat']['variante']!='lu_corrigee'),certificats_rejoues=replays)


if __name__=='__main__':
    print(verifier_archive());print('Sondes :',len(inventaire(SONDES)['fichiers']),'objets vérifiés')

"""Vérification géométrique autonome, schéma public et archives du prototype."""
import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path

if __package__:
    from ._oracle_assemblage import verifier
else:
    from _oracle_assemblage import verifier


SCHEMA = 'vinkulum.assemblage.prototype.1'
SCHEMA_PUBLIC = 'vinkulum.assemblage.1'


def rationnel(s):
    if not isinstance(s, str) or len(s) > 5000:
        raise ValueError('rationnel canonique dans le budget requis')
    q = F(s)
    if str(q) != s or max(q.numerator.bit_length(), q.denominator.bit_length()) > 8192:
        raise ValueError('rationnel non canonique ou hors budget')
    return q


def modele_converti(modele, conversion):
    m = copy.deepcopy(modele)
    for c in m['contraintes']:
        c['pa'] = [conversion(x) for x in c['pa']]
        c['pb'] = [conversion(x) for x in c['pb']]
        if c['type'] == 'distance':
            c['longueur'] = conversion(c['longueur'])
        elif 'ra' in c or 'rb' in c:
            c['ra'] = [list(map(conversion, l)) for l in c['ra']]
            c['rb'] = [list(map(conversion, l)) for l in c['rb']]
    m['jauges'] = [(i, conversion(v)) for i, v in m['jauges']]
    return m


def fabriquer(modele, centre, rayons, inverse):
    if __package__:
        from ._geometrie_certifiee import inclusion
    else:
        from _geometrie_certifiee import inclusion
    inclusion(modele, centre, rayons, inverse)
    enc = lambda x: str(F(x))
    document = {'schema': SCHEMA, 'modele': modele_converti(modele, enc),
                'centre': list(map(enc, centre)), 'rayons': list(map(enc, rayons)),
                'inverse_approche': [list(map(enc, l)) for l in inverse]}
    verifier_document(document)
    return document


def verifier_document(document):
    if isinstance(document, dict) and document.get("schema") == SCHEMA_PUBLIC:
        return verifier_public(document)
    if not isinstance(document, dict) or set(document) != {
            'schema', 'modele', 'centre', 'rayons', 'inverse_approche'} or document['schema'] != SCHEMA:
        raise ValueError('schéma expérimental inconnu')
    try:
        m = modele_converti(document['modele'], rationnel)
        centre = list(map(rationnel, document['centre']))
        rayons = list(map(rationnel, document['rayons']))
        inverse = [list(map(rationnel, l)) for l in document['inverse_approche']]
        return verifier(m, centre, rayons, inverse)
    except (KeyError, TypeError, IndexError, ZeroDivisionError, OverflowError) as e:
        raise ValueError('document géométrique mal formé') from e


def poses_converties(poses, conversion):
    if not isinstance(poses, list) or not 1 <= len(poses) <= 4:
        raise ValueError('nombre de poses hors domaine')
    out = []
    for p in poses:
        if not isinstance(p, dict) or set(p) != {'position', 'position_basse', 'rotation'}:
            raise ValueError('pose mal formée')
        if len(p['position']) != 3 or len(p['position_basse']) != 3 or len(p['rotation']) != 3 or any(len(l) != 3 for l in p['rotation']):
            raise ValueError('dimensions de pose incohérentes')
        out.append({'position': list(map(conversion, p['position'])),
                    'position_basse': list(map(conversion, p['position_basse'])),
                    'rotation': [list(map(conversion, l)) for l in p['rotation']]})
    return out


def fabriquer_public(modele, centre, rayons, inverse, photographie):
    if __package__:
        from ._geometrie_certifiee import bornes_pose
    else:
        from _geometrie_certifiee import bornes_pose
    document = fabriquer(modele, centre, rayons, inverse)
    document['schema'] = SCHEMA_PUBLIC
    document['pose_initiale'] = poses_converties(photographie['poses'], lambda x: str(F(x)))
    bornes = bornes_pose(photographie, centre, rayons)
    document['bornes_pose'] = [{'position': list(map(str, b['position'])),
        'rotation_frobenius_carree': str(b['rotation_frobenius_carree'])} for b in bornes]
    verifier_public(document)
    return document


def verifier_public(document):
    if not isinstance(document, dict) or set(document) != {'schema', 'modele', 'centre',
            'rayons', 'inverse_approche', 'pose_initiale', 'bornes_pose'} or document['schema'] != SCHEMA_PUBLIC:
        raise ValueError('schéma géométrique inconnu')
    if __package__:
        from ._geometrie_certifiee import bornes_pose
    else:
        from _geometrie_certifiee import bornes_pose
    try:
        m = modele_converti(document['modele'], rationnel)
        if any(c['type'] == 'liaison' and not {'ra', 'rb'} <= set(c) for c in m['contraintes']):
            raise ValueError('repères matériels absents du certificat public')
        centre = list(map(rationnel, document['centre']))
        rayons = list(map(rationnel, document['rayons']))
        inverse = [list(map(rationnel, l)) for l in document['inverse_approche']]
        poses = poses_converties(document['pose_initiale'], rationnel)
        if len(poses) != m['corps'] or len(document['bornes_pose']) != len(poses):
            raise ValueError('nombre de corps incohérent')
        bornes = verifier(m, centre, rayons, inverse)
        distances = bornes_pose({'poses': poses}, centre, rayons)
        for annonce, calcule in zip(document['bornes_pose'], distances, strict=True):
            if set(annonce) != {'position', 'rotation_frobenius_carree'} or len(annonce['position']) != 3:
                raise ValueError('bornes de pose mal formées')
            if any(rationnel(a) < c for a, c in zip(annonce['position'], calcule['position'], strict=True)) or rationnel(annonce['rotation_frobenius_carree']) < calcule['rotation_frobenius_carree']:
                raise ValueError('borne de distance à la pose initiale insuffisante')
        return {'dimension': len(centre), 'existence_unique_locale': True,
                'bornes_coordonnees': bornes, 'bornes_pose': distances}
    except (KeyError, TypeError, IndexError, ZeroDivisionError, OverflowError) as e:
        raise ValueError('document géométrique mal formé') from e


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('document', type=Path)
    args = parser.parse_args()
    resultat = verifier_document(json.loads(args.document.read_text()))
    print(json.dumps(resultat, default=str, indent=2))

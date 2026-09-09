"""Taylor validé du pendule : référence rationnelle indépendante du noyau."""
from fractions import Fraction as F
from intervalle_temporel import I, sincos


def etat_de(valeurs):
    if not isinstance(valeurs, (list, tuple)) or len(valeurs) != 2:
        raise ValueError('deux composantes theta, omega requises')
    return list(map(I.de, valeurs))


def ordre_valide(ordre):
    if type(ordre) is not int or not 2 <= ordre <= 12:
        raise ValueError('ordre hors domaine')


def coefficients(etat, k, ordre):
    a, b = [I.de(etat[0])], [I.de(etat[1])]
    s0, c0 = sincos(a[0])
    s, c = [s0], [c0]
    k = I.de(k)
    for n in range(ordre):
        a.append(b[n]/(n+1))
        b.append(-k*s[n]/(n+1))
        s.append(sum(((j+1)*a[j+1]*c[n-j] for j in range(n+1)), I.de(0))/(n+1))
        c.append(-sum(((j+1)*a[j+1]*s[n-j] for j in range(n+1)), I.de(0))/(n+1))
    return a, b


def picard(etat, domaine, h, k):
    h, k = F(h), I.de(k)
    if h <= 0 or h*max(F(1), k.mag()) >= 1:
        raise ValueError('contraction de Picard non établie')
    if any(not domaine[i].lo <= etat[i].lo <= etat[i].hi <= domaine[i].hi for i in range(2)):
        raise ValueError('état initial hors domaine de Picard')
    s, _ = sincos(domaine[0])
    image = [etat[0]+I(0,h)*domaine[1], etat[1]-I(0,h)*k*s]
    if not all(image[i].interieur_de(domaine[i]) for i in range(2)):
        raise ValueError('inclusion stricte de Picard non établie')


def verifier_pas(etat, domaine, h, k, ordre=8):
    ordre_valide(ordre)
    etat, domaine = etat_de(etat), etat_de(domaine)
    h = F(h)
    picard(etat, domaine, h, k)
    polynome = coefficients(etat, k, ordre-1)
    reste = coefficients(domaine, k, ordre)
    return [sum((polynome[i][j]*h**j for j in range(ordre)), I.de(0))
            +reste[i][ordre]*h**ordre for i in range(2)]


def pas(etat, h, k, ordre=8):
    ordre_valide(ordre)
    etat, h, k = etat_de(etat), F(h), I.de(k)
    if h <= 0 or h*max(F(1), k.mag()) >= 1:
        raise ValueError('pas hors domaine')
    s, _ = sincos(etat[0])
    delta = [h*(etat[1].mag()+1), h*((k*s).mag()+1)]
    for _ in range(12):
        domaine = [etat[i]+I(-delta[i],delta[i]) for i in range(2)]
        try:
            resultat = verifier_pas(etat, domaine, h, k, ordre)
            return resultat, domaine
        except ValueError:
            delta = [2*d for d in delta]
    raise ValueError('aucun domaine de Picard établi')

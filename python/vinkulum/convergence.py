"""vinkulum.convergence — « mon pas de temps est-il assez fin ? »

C'est LA question qu'un utilisateur nouveau ne pose pas, et le noyau ne peut
pas y répondre a priori : un schéma dissipatif rend TOUJOURS un résultat, et
un pas mille fois trop grand donne une trajectoire lissée, plausible, et
entièrement fausse. Le critère « h < 1/(20·f_max) » n'aide pas — sur une
structure, f_max vaut 1e5 Hz et le schéma amortit justement ces modes-là :
c'est son rôle.

Ce qui répond, c'est la mesure a posteriori — refaire au pas divisé et
regarder si ça bouge. Tout le monde le fait ; personne ne le fait
systématiquement. Ce module en fait un appel.

Il sert aussi à la question voisine, qui a mordu deux fois sur des bancs
étrangers (6barmech, multibarmech) : **avant de comparer deux solveurs à un
instant donné, mesurer de combien le solveur de RÉFÉRENCE se déplace quand on
raffine SON propre pas.** Si ce déplacement dépasse l'écart entre les deux
codes, la comparaison ne dit rien.
"""

import numpy as np

from ._vinkulum import Noyau  # noqa: F401  (ré-export implicite pour les types)


def etude(monte, t_end, h, facteurs=(1, 2, 4), grandeur=None, rho=0.6):
    """Rejoue la MÊME simulation à h, h/2, h/4… et mesure ce qui bouge.

    `monte` : fonction sans argument rendant un `Noyau` neuf, déjà assemblé.
    `grandeur` : fonction(Noyau) -> vecteur ; par défaut, toutes les positions.

    Rend une liste de (h, écart au pas le plus fin, ordre observé).
    Les facteurs doivent être positifs et strictement croissants (au moins
    deux). L'ordre n'est estimé que pour trois pas successifs de même rapport,
    avec des écarts non nuls ; sinon il vaut None. Les observables sont copiées.
    Un modèle indépendant au même instant initial est requis à chaque appel.
    Cette étude mesure le raffinement de l'observable finale : elle ne prouve
    ni une erreur absolue, ni la justesse physique du modèle.
    """
    t_end, h, rho = float(t_end), float(h), float(rho)
    if not np.isfinite(t_end) or not np.isfinite(h) or h <= 0:
        raise ValueError("date finale finie et pas fini > 0 requis")
    if not np.isfinite(rho):
        raise ValueError("rho doit être fini")
    facteurs = np.array(tuple(facteurs), dtype=float)
    if (facteurs.ndim != 1 or facteurs.size < 2
            or not np.isfinite(facteurs).all() or (facteurs <= 0).any()
            or not (facteurs[1:] > facteurs[:-1]).all()):
        raise ValueError("au moins deux facteurs finis, positifs et strictement croissants requis")
    with np.errstate(over='ignore', under='ignore'):
        pas = h / facteurs
    if (not np.isfinite(pas).all() or (pas <= 0).any()
            or not (pas[1:] < pas[:-1]).all()):
        raise ValueError("pas raffinés non représentables")
    if grandeur is None:
        def grandeur(n):
            return np.asarray(n.etat()[1]).ravel()
    vals = []
    modeles = []  # garder les références empêche aussi la réutilisation des id
    t0 = None
    for hh in pas:
        n = monte()
        if not isinstance(n, Noyau):
            raise ValueError("monte doit rendre un Noyau")
        if any(n is precedent for precedent in modeles):
            raise ValueError("monte doit rendre un modèle indépendant à chaque appel")
        debut = n.t()
        if not np.isfinite(debut) or t_end <= debut:
            raise ValueError("la date finale doit suivre l'instant initial fini")
        if t0 is not None and debut != t0:
            raise ValueError("les modèles doivent avoir le même instant initial")
        t0 = debut
        modeles.append(n)
        n.simule(t_end, float(hh), rho=rho, tous=10 ** 9)
        if n.t() != t_end:
            raise RuntimeError("simulation incomplète : date finale demandée non atteinte")
        v = np.array(grandeur(n), dtype=float, copy=True)
        if not v.size or not np.isfinite(v).all():
            raise ValueError("observable non vide et finie requise")
        if vals and v.shape != vals[0][1].shape:
            raise ValueError("dimensions des observables incohérentes")
        vals.append((float(hh), v))
    ref = vals[-1][1]
    # L'ORDRE SE MESURE SUR DES DIFFERENCES SUCCESSIVES, PAS SUR L'ECART A LA
    # REFERENCE. Le pas le plus fin n'est pas la solution exacte : quand
    # l'ecart s'en approche, il est domine par l'erreur de la reference
    # elle-meme et l'ordre apparent part n'importe ou (mesure : 5,21 puis 1,48
    # puis 1,76 sur un cas d'ordre 2 franc). Richardson compare v(h) a v(h/r)
    # et ne suppose aucune solution exacte.
    def ecart(a, b):
        with np.errstate(over='ignore', invalid='ignore'):
            e = float(np.max(np.abs(a - b)))
        if not np.isfinite(e):
            raise ValueError("écart entre observables non fini")
        return e

    d = [ecart(vals[i][1], vals[i + 1][1]) for i in range(len(vals) - 1)]
    out = []
    for i, (hh, v) in enumerate(vals[:-1]):
        e = ecart(v, ref)
        ordre = None
        if i + 1 < len(d) and d[i + 1] > 0 and d[i] > 0:
            log_r1 = np.log(facteurs[i + 1]) - np.log(facteurs[i])
            log_r2 = np.log(facteurs[i + 2]) - np.log(facteurs[i + 1])
            if log_r1 > 0 and np.isclose(log_r1, log_r2, rtol=1e-12, atol=0):
                estimation = float((np.log(d[i]) - np.log(d[i + 1])) / log_r1)
                if np.isfinite(estimation):
                    ordre = estimation
        out.append((hh, e, ordre))
    return out


def verdict(monte, t_end, h, tol, **kw):
    """Accord du calcul au premier pas avec le plus fin à `tol` près.

    Rend (verdict, écart mesuré, table), pour l'observable finale choisie.
    tol est absolue, finie et >=0, dans les unités de cette observable.
    Ce verdict ne garantit ni une erreur absolue ni la validité du modèle.
    """
    tol = float(tol)
    if not np.isfinite(tol) or tol < 0:
        raise ValueError("tol doit être finie et >= 0")
    t = etude(monte, t_end, h, **kw)
    return (t[0][1] <= tol), t[0][1], t


def demo():
    """Autotest, en DEUX cas — et leur différence est l'enseignement.

    (a) **une seule échelle** : pendule rigide. L'ordre observé doit valoir 2,
        celui du schéma. C'est ce qui valide la MESURE d'ordre elle-même.

    (b) **deux échelles séparées** : poutre raide (flexion à 17 Hz, traction à
        500 Hz). À pas grossier le mode rapide n'est pas résolu du tout — le
        schéma l'amortit, ce qui est son rôle — et l'ordre apparent n'est
        PLUS 2. Ce n'est pas un défaut : c'est l'information. Un ordre observé
        loin de 2 dit « ce pas ne résout pas tout ce que le modèle contient ».
    """
    from ._vinkulum import Noyau

    J = [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]

    def pendule():
        n = Noyau([0.0, 0.0, -9.81])
        c = n.corps("c", 1.0, J, [0.5, 0, 0])
        n.liaison("p", None, c, pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 2])
        return n

    def raide():
        n = Noyau([0.0, 0.0, -9.81])
        a = n.corps("a", 1.0, J, [0, 0, 0])
        b = n.corps("b", 1.0, J, [1, 0, 0])
        n.poutre("p", a, b, 1e7, 5e6, 2e3, 4e3)
        n.liaison("l", None, a, pa=[0, 0, 0])
        return n

    print("╔═ vinkulum — convergence : « mon pas est-il assez fin ? »")
    ta = etude(pendule, 1.0, 8e-3, facteurs=(1, 2, 4, 8))
    print("║ (a) une seule echelle — pendule rigide")
    for hh, e, o in ta:
        print(f"║     h = {hh:.1e} : ecart {e:.3e} m" + (f"   ordre {o:.2f}" if o else ""))
    tb = etude(raide, 0.2, 4e-3, facteurs=(1, 2, 4, 8, 16))
    print("║ (b) deux echelles separees — poutre raide (flexion 17 Hz, traction 500 Hz)")
    for hh, e, o in tb:
        print(f"║     h = {hh:.1e} : ecart {e:.3e} m" + (f"   ordre {o:.2f}" if o else ""))
    ok, e, _ = verdict(raide, 0.2, 4e-3, tol=1e-6)
    ok2, e2, _ = verdict(raide, 0.2, 5e-5, tol=1e-6)

    # ET LE NOYAU SAIT DEJA LE FAIRE — il ne le fait simplement pas de
    # lui-meme. Le pas ADAPTATIF (resid7u de demi-pas, Hibbitt & Karlsson 1979)
    # part du pas qui ment et descend jusqu'a ce que la tolerance soit tenue.
    # On ne change pas le defaut : ce serait changer la reproductibilite de
    # tout modele existant, et le pas est un choix de l'utilisateur. Mais il
    # doit SAVOIR que ca existe et ce que ca vaut, et ca se mesure.
    ref = raide()
    ref.simule(0.2, 5e-5, tous=10 ** 9)
    zref = ref.etat()[1][-1][2]
    ada = []
    for tolr in (1e-3, 1e-5):
        n = raide()
        n.simule(0.2, 4e-3, adaptatif=tolr, tous=10 ** 9)
        ada.append((tolr, abs(n.etat()[1][-1][2] - zref), n.adapt_stats()[1]))
    nf = raide()
    nf.simule(0.2, 4e-3, tous=10 ** 9)
    e_fixe = abs(nf.etat()[1][-1][2] - zref)
    print(f"║ le meme depart h = 4e-3, en ADAPTATIF : "
          + " · ".join(f"tol {t:.0e} → ecart {d:.1e} (pas final {p:.1e})"
                       for t, d, p in ada))
    print(f"║ contre {e_fixe:.1e} a pas fixe : le noyau SAIT corriger, il ne le fait")
    print("║ pas de lui-meme — le defaut n'est pas touche, mais il est mesure.")
    assert ada[-1][1] < e_fixe / 100, ("l'adaptatif n'ameliore pas", ada, e_fixe)
    print(f"║ verdict a h = 4e-3 (tol 1e-6 m) : {'ASSEZ FIN' if ok else 'TROP GROSSIER'}"
          f" (ecart {e:.2e})   ·   a h = 5e-5 : {'ASSEZ FIN' if ok2 else 'TROP GROSSIER'}"
          f" (ecart {e2:.2e})")

    # · (a) valide la MESURE : sur une seule echelle, l'ordre observe est celui
    #   du schema. Sans ce controle, un ordre bizarre en (b) serait ambigu.
    oa = [o for _, _, o in ta if o is not None]
    assert 1.7 < np.median(oa) < 2.3, ("l'ordre du schema n'est pas retrouve", oa)
    # · (b) le raffinement rapproche, et l'ordre apparent N'EST PAS 2 : le pas
    #   grossier ne resout pas le mode rapide, et la table le dit.
    assert tb[0][1] > 100 * tb[-1][1], ("le raffinement ne rapproche pas", tb)
    ob = [o for _, _, o in tb if o is not None]
    assert max(ob) > 3.0, ("l'ordre apparent devrait sortir de 2 sur un cas raide", ob)
    # · le verdict discrimine
    assert not ok and ok2, ("le verdict ne discrimine pas", e, e2)
    print("╚ convergence OK — un ordre observe loin de 2 dit « ce pas ne resout pas tout »")
    return dict(pendule=ta, raide=tb, grossier=e, fin=e2)


if __name__ == "__main__":
    demo()

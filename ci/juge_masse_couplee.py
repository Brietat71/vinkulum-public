"""Juge de champs couplés : facteur massique audité, six charges et opérateurs.

Le facteur du juge est indépendant de la réduction. Ses erreurs de modèle
sont encadrées en rationnels ; la qualification reste numérique, comme les
QR et SVD du juge historique. Aucun Gram de petite combinaison de charge.
"""
CIBLE=1e-6


def juger_couple(d, m, metrique, candidat, reference, *, bibliotheque):
    """Erreurs physiques par charge et sur toutes leurs combinaisons réelles.

    Les six colonnes de référence sont normalisées avant QR pour éviter
    qu'une différence d'amplitude entre axial et flexion soit prise pour
    une perte de rang. La même transformation s'applique à l'erreur :
    ||B_erreur R^-1|| reste la norme opérateur relative exacte. Le rang
    numérique insuffisant est refusé ; aucune direction n'est supprimée.
    """
    import numpy as np
    from scipy.linalg import cholesky, qr, solve_triangular, svdvals
    from scipy.sparse import issparse
    if any(np.iscomplexobj(a) for a in (d, m, metrique, candidat, reference)):
        raise ValueError("le juge physique requiert des données réelles")
    candidat, reference = np.asarray(candidat), np.asarray(reference)
    metrique = np.asarray(metrique, dtype=float)
    if (candidat.ndim != 3 or candidat.shape != reference.shape
            or candidat.shape[0] == 0 or candidat.shape[1] < 6 or candidat.shape[2] != 6):
        raise ValueError("champs de mêmes dimensions fréquences×coordonnées×6 requis")
    total = candidat.shape[1]
    if (len(d.shape) != 2 or d.shape[1] != total or m.shape != (total, total)
            or metrique.shape != (6, 6)):
        raise ValueError("dimensions incompatibles de D, M ou de la métrique de port")
    for nom, a in (("D", d), ("M", m), ("métrique", metrique),
                   ("candidat", candidat), ("référence", reference)):
        valeurs = a.data if issparse(a) else np.asarray(a)
        if not np.all(np.isfinite(valeurs)):
            raise ValueError(nom+" contient des valeurs non finies")
    from racine_masse_juge import RacineMasseJuge
    racine = RacineMasseJuge(m,bibliotheque=bibliotheque)
    eps = np.finfo(float).eps
    if np.max(np.abs(metrique-metrique.T)) > 64*eps*np.max(np.abs(metrique)):
        raise ValueError("métrique de port non symétrique")
    # Même convention que la normalisation des ports : Cholesky du
    # triangle inférieur ; une dissymétrie significative est refusée.
    l = cholesky(metrique, lower=True)
    erreurs = {n: [] for n in ("masse", "deformation", "port")}
    operateurs = {n: [] for n in erreurs}
    for frequence, (x, y) in enumerate(zip(candidat, reference, strict=True)):
        delta = x-y
        couples = (("masse", racine.r@delta, racine.r@y),
                   ("deformation", d@delta, d@y),
                   ("port", l.T@delta[-6:], l.T@y[-6:]))
        for nom, be, br in couples:
            if not np.all(np.isfinite(be)) or not np.all(np.isfinite(br)):
                raise ArithmeticError(f"champ pondéré non fini ({nom}, fréquence {frequence})")
            numerateur = np.linalg.norm(be, axis=0)
            denominateur = np.linalg.norm(br, axis=0)
            if (not np.all(np.isfinite(numerateur)) or not np.all(np.isfinite(denominateur))
                    or np.any(denominateur <= 0)):
                raise ArithmeticError(f"référence de norme nulle ou non représentable ({nom}, fréquence {frequence})")
            colonnes = numerateur/denominateur
            reference_normalisee = br/denominateur[None, :]
            erreur_normalisee = be/denominateur[None, :]
            if br.shape[0] < 6:
                raise ArithmeticError(f"référence de rang inférieur à six ({nom}, fréquence {frequence})")
            _, r = qr(reference_normalisee, mode="economic", check_finite=False)
            valeurs = svdvals(r, check_finite=False)
            seuil_rang = eps*max(reference_normalisee.shape)*valeurs[0]
            if valeurs[-1] <= seuil_rang:
                raise ArithmeticError(f"référence numériquement de rang inférieur à six ({nom}, fréquence {frequence})")
            # Résolution à droite ; aucune inverse explicite ni Gram formé.
            relatif = solve_triangular(r.T, erreur_normalisee.T, lower=True,
                                       check_finite=False).T
            operateur = float(svdvals(relatif, check_finite=False)[0])
            if not np.all(np.isfinite(colonnes)) or not np.isfinite(operateur):
                raise ArithmeticError(f"erreur relative non finie ({nom}, fréquence {frequence})")
            if nom == 'masse':
                colonnes = np.where(colonnes>0,np.nextafter(colonnes*racine.correction,np.inf),0.)
                operateur = float(np.nextafter(operateur*racine.correction,np.inf)) if operateur>0 else 0.
            erreurs[nom].append(colonnes.tolist())
            operateurs[nom].append(operateur)
    maxima = {n: float(np.max(v)) for n, v in erreurs.items()}
    maxima_operateurs = {n: float(np.max(v)) for n, v in operateurs.items()}
    return dict(audit_masse=racine.audit, erreurs=erreurs, maxima=maxima, operateurs=operateurs,
                maxima_operateurs=maxima_operateurs,
                accepte=all(v <= CIBLE for v in (*maxima.values(), *maxima_operateurs.values())))



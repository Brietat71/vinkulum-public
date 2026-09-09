"""vinkulum.urdf — lire un modèle qu'on n'a PAS écrit.

Le verrou du « noyau généraliste » n'est pas la validation : dix modèles
étrangers ont été repris sans trouver une seule erreur de physique. C'est
qu'il n'existait AUCUN moyen d'entrer un modèle autrement qu'en écrivant du
Python contre l'API — donc « poser son modèle dedans » voulait dire le
RÉÉCRIRE, et à ce moment-là ce n'est plus son modèle.

URDF est le format de la robotique (ROS, Bullet, Gazebo, MuJoCo par
conversion). Des milliers de modèles publics existent. Ce module en fait des
`Noyau` sans que rien ne soit retranscrit à la main.

CE QU'IL COUVRE : `link` avec `inertial`, `joint` de type `fixed`,
`revolute`, `continuous`, `prismatic`, `floating`, `planar`, la convention
`rpy` de l'URDF (R = Rz·Ry·Rx), l'inertie exprimée dans son propre repère.

CE QU'IL NE COUVRE PAS, déclaré : `visual` et `collision` (aucune géométrie
n'entre — le noyau n'en a pas besoin pour la dynamique, et les mailles
référencées sont des fichiers externes) ; `mimic` ; les `limit` (le noyau n'a
pas de butée d'articulation générique ; `LoiCouple::Butee` existe mais se pose
à la main) ; `transmission` et `gazebo`. Un `xacro` non déplié est refusé.
"""

import xml.etree.ElementTree as ET

import numpy as np

from ._vinkulum import Noyau


def _fini(v, quoi):
    """Refuse NaN et inf. Un fichier tiers en contient parfois, et un modèle
    à `nan` kg s'intègre SANS ERREUR en produisant du `nan` partout — le pire
    silence possible (mesuré le 3 sept. en corrompant des modèles publics)."""
    if not np.all(np.isfinite(np.asarray(v, float))):
        raise ValueError(f"{quoi} : valeur non finie ({v})")
    return v


def _rpy(s):
    """R = Rz(yaw)·Ry(pitch)·Rx(roll) — la convention URDF, pas une autre."""
    r, p, y = (float(x) for x in (s or "0 0 0").split())
    cr, sr, cp, sp, cy, sy = (np.cos(r), np.sin(r), np.cos(p),
                              np.sin(p), np.cos(y), np.sin(y))
    return np.array([[cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
                     [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
                     [-sp, cp * sr, cp * cr]])


def _origine(e):
    if e is None:
        return np.zeros(3), np.eye(3)
    o = e.find("origin")
    if o is None:
        return np.zeros(3), np.eye(3)
    xyz = np.array([float(x) for x in (o.get("xyz") or "0 0 0").split()])
    return xyz, _rpy(o.get("rpy"))


def _inertie(e):
    """(masse, CdM local, R local du repère d'inertie, tenseur dans CE repère)."""
    ine = e.find("inertial")
    if ine is None:
        return 0.0, np.zeros(3), np.eye(3), np.zeros((3, 3))
    m = _fini(float(ine.find("mass").get("value")), "masse")
    c, R = _origine(ine)
    i = ine.find("inertia")
    v = {k: float(i.get(k, 0.0)) for k in
         ("ixx", "ixy", "ixz", "iyy", "iyz", "izz")}
    J = np.array([[v["ixx"], v["ixy"], v["ixz"]],
                  [v["ixy"], v["iyy"], v["iyz"]],
                  [v["ixz"], v["iyz"], v["izz"]]])
    _fini(J.ravel(), "tenseur d'inertie")
    _fini(c, "centre de masse")
    return m, c, R, J


# type URDF -> (composantes de translation bloquées, de rotation bloquées),
# « axe » = l'indice de l'axe du joint dans son propre repère, mis en 0
_BLOQUE = {"fixed": ([0, 1, 2], [0, 1, 2]),
           "revolute": ([0, 1, 2], [1, 2]),
           "continuous": ([0, 1, 2], [1, 2]),
           "prismatic": ([1, 2], [0, 1, 2]),
           "floating": ([], []),
           "planar": ([0], [1, 2]),
           # `spherical` n'est pas dans la spécification URDF, mais Bullet
           # l'écrit (son `humanoid.urdf`) et c'est une ROTULE — que le noyau
           # fait depuis toujours. Le refuser était un refus injustifié, le
           # même que `ball` côté MJCF.
           "spherical": ([0, 1, 2], [])}


def _repere_axe(a):
    """Repère dont le PREMIER axe est `a` — les liaisons du noyau bloquent des
    indices, pas des directions ; on tourne le repère de liaison à la place."""
    e1 = np.asarray(a, float)
    n = np.linalg.norm(e1)
    if n < 1e-12:
        e1 = np.array([1.0, 0.0, 0.0])
    else:
        e1 = e1 / n
    t = np.array([0.0, 0.0, 1.0]) if abs(e1[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e2 = np.cross(e1, t)
    e2 /= np.linalg.norm(e2)
    return np.column_stack([e1, e2, np.cross(e1, e2)])


def charge(chemin, g=(0.0, 0.0, -9.81), encastre_racine=True, m_min=1e-9,
           butees=True, garde_deg=1.0):
    """Construit un `Noyau` depuis un fichier URDF.

    Rend (noyau, infos) où `infos` porte les correspondances nom→indice, la
    masse totale, et la liste des simplifications APPLIQUÉES — un import qui
    tait ce qu'il a dû arranger est un import qui ment.

    `butees` : pose les `<limit>` des joints `revolute` comme butées de couple.
    La RAIDEUR est dérivée du fichier, pas posée — `k = effort / garde`, où
    `garde_deg` est le dépassement (en degrés) auquel le couple maximal admis
    par l'URDF est atteint. C'est un choix DÉCLARÉ : l'URDF donne un couple
    limite, pas une raideur de butée, et il faut bien la fabriquer.
    Les `prismatic` sont laissés SANS butée — le noyau n'a pas de butée en
    force, seulement en couple — et l'import le dit dans ses notes.
    """
    src = open(chemin, "rb").read()
    # ⚠ `xmlns:xacro` en attribut de <robot> NE veut PAS dire non déplié : les
    # fichiers générés par xacro le gardent presque toujours. Ce qui signe un
    # xacro brut, c'est une substitution `${…}` ou un ÉLÉMENT `<xacro:…>`.
    # Première version : crié au loup sur le premier vrai fichier.
    if b"${" in src or b"<xacro:" in src:
        raise ValueError(f"{chemin} : xacro non déplié — lancer `xacro` d'abord")
    # substitutions ROS non résolues : `$(find …)`, `$(optenv …)`, `$(arg …)`.
    # Le fichier n'est pas un URDF, c'est un gabarit — le dire, plutôt que
    # d'échouer plus loin sur « could not convert string to float ».
    if b"$(" in src:
        raise ValueError(f"{chemin} : substitution ROS non résolue ($(find/optenv/arg)) "
                         "— ce fichier est un gabarit, pas un URDF")
    # ⚠ DES OCTETS NULS EN FIN DE FICHIER. Le `humanoid.urdf` de Bullet en
    # porte un après `</robot>` : XML invalide au sens strict, et ElementTree
    # rend « not well-formed (invalid token) » sans dire où le mal est. C'est
    # un artefact d'écriture fréquent, pas une erreur de modèle — on le
    # nettoie et on le DIT, plutôt que de refuser un fichier utilisable.
    nuls = src.count(b"\x00")
    if nuls:
        src = src.replace(b"\x00", b"")
    try:
        root = ET.fromstring(src)
    except ET.ParseError as ex:
        raise ValueError(f"{chemin} : XML invalide — {ex}") from ex
    if root.tag != "robot":
        raise ValueError(f"{chemin} : racine <{root.tag}>, attendu <robot>")

    notes = []
    links = {e.get("name"): e for e in root.findall("link")}
    joints = [e for e in root.findall("joint")]
    # SENTINELLE D'OUBLI : on ne lit que les enfants DIRECTS de <robot>. Un
    # <joint> ou un <link> imbriqué ailleurs serait ignoré en SILENCE, et un
    # import incomplet est pire qu'un import refusé.
    # ⚠ ON NE COMPTE QUE LES DÉFINITIONS. Un `<transmission>` contient des
    # `<joint name=...>` qui REFERENCENT une articulation pour lui attacher
    # une interface matérielle ; ce ne sont pas des définitions, et les
    # compter faisait refuser des modèles publics parfaitement valides
    # (husky et racecar de Bullet, mesuré le 3 sept.). Une définition de
    # joint a un <parent> et un <child> ; une définition de link a un nom.
    def _defs(tag):
        return [e for e in root.iter(tag)
                if (tag == "link" and e.get("name"))
                or (tag == "joint" and e.find("parent") is not None
                    and e.find("child") is not None)]

    for tag, pris in (("joint", joints), ("link", list(links.values()))):
        tous = _defs(tag)
        if len(tous) != len(pris):
            raise ValueError(
                f"{chemin} : {len(tous)} <{tag}> definis dans le document mais "
                f"{len(pris)} enfants directs de <robot> — le reste serait ignoré")
    parent = {}
    for j in joints:
        c = j.find("child").get("link")
        if c in parent:
            raise ValueError(f"link « {c} » a deux parents — URDF exige un ARBRE")
        parent[c] = (j.find("parent").get("link"), j)

    if nuls:
        notes.append(f"{nuls} octet(s) NUL retire(s) du fichier — XML invalide "
                     "au sens strict, artefact d'ecriture")
    racines = [n for n in links if n not in parent]
    if len(racines) != 1:
        raise ValueError(f"{chemin} : {len(racines)} racines ({racines}) — attendu une")

    # poses des repères de link, propagées depuis la racine
    pose = {racines[0]: (np.zeros(3), np.eye(3))}
    reste = [j for j in joints]
    while reste:
        avance = False
        for j in list(reste):
            p, c = j.find("parent").get("link"), j.find("child").get("link")
            if p in pose:
                tp, Rp = pose[p]
                to, Ro = _origine(j)
                pose[c] = (tp + Rp @ to, Rp @ Ro)
                reste.remove(j)
                avance = True
        if not avance:
            # DEUX causes bloquent la propagation, et le message doit dire
            # LAQUELLE : un parent qui n'existe pas dans <link>, ou un vrai
            # cycle. Accuser le cycle dans les deux cas envoie chercher au
            # mauvais endroit (sorti par le banc lecteurs_casses).
            orphelins = sorted({j.find("parent").get("link") for j in reste}
                               - set(links) - set(pose))
            if orphelins:
                raise ValueError(
                    f"{chemin} : parent(s) inexistant(s) dans <link> : "
                    f"{orphelins} — référencés par un <joint>")
            raise ValueError(
                f"{chemin} : cycle dans l'arbre des links — "
                f"{len(reste)} joint(s) non propageable(s)")

    n = Noyau(list(g))
    idx, mtot = {}, 0.0
    for nom, e in links.items():
        m, c, Ri, J = _inertie(e)
        tl, Rl = pose[nom]
        if m <= 0.0:
            # licite en URDF (repère intermédiaire) et impossible en multicorps :
            # un corps sans masse rend la matrice singulière.
            m, J = m_min, np.eye(3) * m_min
            Ri = np.eye(3)
            notes.append(f"{nom} : sans inertie dans l'URDF — masse {m_min:g} kg imposée")
        elif np.min(np.linalg.eigvalsh(J)) <= 0.0:
            # ⚠ INERTIE NULLE AVEC UNE MASSE NON NULLE. Le laikago de Bullet
            # écrit `ixx=iyy=izz=0` pour un châssis de 13,7 kg : le fichier est
            # physiquement incomplet, et rien dans l'URDF ne permet de la
            # reconstruire (aucune géométrie n'est lue). On impose le minimum
            # qui rende le modèle intégrable et on le DIT — la dynamique de
            # ROTATION de ce corps n'est alors pas celle du fichier.
            J = np.eye(3) * max(m * 1e-6, m_min)
            Ri = np.eye(3)
            notes.append(f"{nom} : <inertia> nulle pour {m:.3f} kg — inertie minimale "
                         f"imposée, la dynamique de ROTATION de ce corps n'est PAS "
                         f"celle du fichier")
        Rc = Rl @ Ri
        idx[nom] = n.corps(nom, m, list(J.flatten()), list(tl + Rl @ c),
                           rot=list(Rc.flatten()))
        mtot += m

    for j in joints:
        typ = j.get("type")
        if typ not in _BLOQUE:
            raise ValueError(f"joint « {j.get('name')} » : type « {typ} » non couvert")
        bt, br = _BLOQUE[typ]
        p, c = j.find("parent").get("link"), j.find("child").get("link")
        ax = j.find("axis")
        a = [float(x) for x in (ax.get("xyz") if ax is not None else "1 0 0").split()]
        _fini(a, f"axe du joint « {j.get('name')} »")
        if typ in ("revolute", "continuous", "prismatic") and np.linalg.norm(a) < 1e-12:
            # un axe NUL ne définit aucune liaison : le degré libéré serait
            # arbitraire. Refuser, plutôt que d'en choisir un au hasard.
            raise ValueError(f"joint « {j.get('name')} » de type {typ} : axe NUL")
        tj, Rj = pose[c]                       # le repère du joint = celui de l'enfant
        Ra = Rj @ _repere_axe(a)               # premier axe = l'axe du joint
        # tout est exprimé dans le repère du CORPS parent, pas du link parent
        tp, Rp = pose[p]
        mp, cp, Rip, _ = _inertie(links[p])
        rp = tp + Rp @ (cp if mp > 0 else np.zeros(3))
        Rcp = Rp @ (Rip if mp > 0 else np.eye(3))
        n.liaison(j.get("name"), idx[p], idx[c],
                  pa=list(Rcp.T @ (tj - rp)), ra=list((Rcp.T @ Ra).flatten()),
                  bloque_t=bt, bloque_r=br)

    if butees:
        for j in joints:
            lim = j.find("limit")
            typ = j.get("type")
            if lim is None or typ not in ("revolute", "prismatic"):
                continue
            lo, hi = float(lim.get("lower", 0.0)), float(lim.get("upper", 0.0))
            if not hi > lo:
                continue
            if typ == "prismatic":
                notes.append(f"{j.get('name')} : <limit> ignore — le noyau n'a pas "
                             "de butee en FORCE, seulement en couple")
                continue
            eff = float(lim.get("effort", 0.0))
            if eff <= 0.0:
                notes.append(f"{j.get('name')} : <limit> sans effort — butee non posee")
                continue
            k = eff / np.radians(garde_deg)
            ax = j.find("axis")
            a = [float(x) for x in (ax.get("xyz") if ax is not None else "1 0 0").split()]
            p, c = j.find("parent").get("link"), j.find("child").get("link")
            tj, Rj = pose[c]
            _, Rp = pose[p]
            mp, _, Rip, _ = _inertie(links[p])
            axe_a = (Rp @ (Rip if mp > 0 else np.eye(3))).T @ (Rj @ np.asarray(a, float))
            n.couple(f"{j.get('name')}_butee", idx[p], idx[c], list(axe_a),
                     ("butee", [k, 0.02 * k, lo, hi]))

    if encastre_racine:
        r = racines[0]
        tl, Rl = pose[r]
        n.liaison("_bati", None, idx[r], pa=list(tl), ra=list(Rl.flatten()))

    # la masse importée doit être celle du fichier, aux corps sans inertie près
    m_fichier = sum(_inertie(e)[0] for e in links.values())
    if abs(mtot - m_fichier) > 1e-6 * max(m_fichier, 1.0) + m_min * len(links) * 1.01:
        raise ValueError(f"{chemin} : masse importée {mtot} ≠ masse du fichier {m_fichier}")
    return n, dict(noms=idx, racine=racines[0], masse=mtot, notes=notes,
                   butees=sum(1 for j in joints
                              if butees and j.find("limit") is not None
                              and j.get("type") == "revolute"
                              and float(j.find("limit").get("effort", 0.0)) > 0.0),
                   m_fichier=m_fichier, n_links=len(links), n_joints=len(joints),
                   n_corps=len(links), n_liaisons=len(joints),
                   types=sorted({j.get("type") for j in joints}))


def _conservation(n, graine=7, t=0.12, h=1e-3):
    """dL/L sur un système LIBRE et SANS gravité — un contrôle EXACT.

    Sans effort extérieur, le moment cinétique est conservé quelles que soient
    les liaisons internes. C'est la seule propriété qu'on puisse exiger d'un
    modèle dont on ne connaît pas la solution, et elle vaut donc pour
    n'importe quel fichier qu'un tiers apporterait.

    L'ÉNERGIE ne peut pas servir de critère strict ici : le schéma
    α-généralisé DISSIPE volontairement les hautes fréquences (ρ∞ < 1), et sur
    ces modèles la dérive atteint 6e-4 sans que rien ne soit faux.
    """
    import numpy as np
    st = n.etat()
    rng = np.random.default_rng(graine)
    v = [list(rng.normal(0, 0.3, 3)) for _ in st[1]]
    w = [list(rng.normal(0, 0.5, 3)) for _ in st[1]]
    n.pose_etat(0.0, st[1], st[2], v, w)
    try:
        n.assemble()
    except Exception:                      # noqa: BLE001 — pose déjà admissible
        pass
    l0 = np.asarray(n.moment())
    n.simule(t, h, tous=10 ** 9)
    l1 = np.asarray(n.moment())
    return float(np.linalg.norm(l1 - l0) / max(np.linalg.norm(l0), 1e-12))

def demo():
    """Autotest sur DEUX robots publics, téléchargés tels quels — le KUKA iiwa
    et le R2D2 de Bullet. Rien n'y est de moi que le lecteur.
    """
    import os
    d = os.path.join(os.path.dirname(__file__), "donnees")
    print("╔═ vinkulum — URDF : lire un modèle qu'on n'a pas écrit")
    out = {}
    for f, attendu in (("kuka_iiwa.urdf", 8), ("r2d2.urdf", 16)):
        n, info = charge(os.path.join(d, f))
        phi = max(abs(x) for x in n.phi())
        n.simule(0.2, 1e-3, tous=10 ** 9)
        phi2 = max(abs(x) for x in n.phi())
        e = n.etat()[1]
        print(f"║ {f:16s} {info['n_links']:2d} corps · {info['n_joints']:2d} liaisons · "
              f"{info['masse']:7.3f} kg · |Φ| {phi:.0e} → {phi2:.0e} après 200 pas")
        print(f"║     types de liaison : {', '.join(info['types'])} · "
              f"{info['butees']} butee(s) d'articulation posee(s)")
        for s in info["notes"][:2]:
            print(f"║     simplification : {s}")
        assert info["n_links"] == attendu, (f, info["n_links"])
        assert phi < 1e-12 and phi2 < 1e-9, ("contraintes", f, phi, phi2)
        assert all(np.all(np.isfinite(p)) for p in e), ("etat non fini", f)
        out[f] = info

    # LES BUTEES DOIVENT MORDRE. Les poser ne prouve rien : on lache le KUKA
    # sous gravite avec et sans, et la trajectoire doit differer. Sans ce
    # controle, un import qui les ignorerait passerait le banc.
    poses = {}
    # h = 1e-3, MESURÉ (7 sept.) : à 2e-3 cette chute avec butées en pénalité est
    # un tirage à pile ou face — elle passait ou lâchait selon le chemin de Newton
    # (ordre 0 ou 1, règles près du seuil), et à 1,8 / 1,9 / 2,1 / 2,2e-3 elle
    # lâchait aussi avec le solveur d'origine. De 1,0 à 1,9e-3 elle passe dans
    # tous les modes essayés. La chute est chaotique (écart avec/sans butées de
    # 0,44 à 0,68 m selon le chemin) : ce banc ne juge que le SIGNE de l'écart.
    for b in (False, True):
        n, _ = charge(os.path.join(d, "kuka_iiwa.urdf"), butees=b)
        n.simule(3.0, 1e-3, tous=10 ** 9)
        poses[b] = np.array(n.etat()[1])
    ecart = float(np.max(np.abs(poses[True] - poses[False])))
    print(f"║ butees : le KUKA lache 3 s sous gravite s'ecarte de {ecart:.3f} m "
          f"selon qu'elles sont posees ou non")
    assert ecart > 0.05, ("les butees d'articulation ne mordent pas", ecart)
    out["butees_ecart"] = ecart
    # ── LE LOT, comme pour MJCF : un taux, pas des cas choisis ──
    du = os.path.join(d, "urdfs")
    lot, refus = [], []
    for f in sorted(os.listdir(du)):
        if not f.endswith(".urdf"):
            continue
        try:
            nl, il = charge(os.path.join(du, f))
            nl.simule(0.05, 1e-3, tous=10 ** 9)
            assert all(np.all(np.isfinite(p)) for p in nl.etat()[1])
            nl2, _ = charge(os.path.join(du, f), g=(0, 0, 0), encastre_racine=False)
            dl = _conservation(nl2)
            lot.append((f, il["n_links"], il["masse"], len(il["notes"]), dl))
        except ValueError as ex:                     # refus MOTIVÉ du lecteur
            refus.append((f, str(ex)[:60]))
    print(f"║ ── le LOT : {len(lot)} robots publics charges, {len(refus)} refuses avec motif ──")
    print("║   " + " · ".join(f"{f[:-5]} {c}c {m:.1f}kg" + (f" ({k}n)" if k else "")
                          for f, c, m, k, _ in lot))
    print(f"║   conservation du moment cinetique (libre, sans gravite) : dL/L max "
          f"{max(x[4] for x in lot):.1e} sur les {len(lot)} modeles")
    assert max(x[4] for x in lot) < 1e-6, ("un modele public ne conserve pas L", lot)
    for f, e in refus:
        print(f"║   refus {f} : {e}")
    assert len(lot) >= 8, ("le lot URDF a maigri", lot, refus)
    out["lot"] = lot

    print("╚ URDF OK — robots publics chargés et intégrés sans une ligne retranscrite")
    return out


if __name__ == "__main__":
    demo()

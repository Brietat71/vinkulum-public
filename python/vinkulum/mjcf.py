"""vinkulum.mjcf — le SECOND format d'entrée, et il ne ressemble pas au premier.

L'URDF décrit un arbre par RÉFÉRENCE (chaque joint nomme son parent et son
enfant) et donne les inerties. MJCF (MuJoCo) décrit un arbre par
IMBRICATION, autorise PLUSIEURS articulations sur un même corps, hérite des
attributs par CLASSES, et surtout **ne donne pas les masses** : elles se
calculent depuis la géométrie et une masse volumique.

Lire deux formats de sémantiques différentes est ce qui distingue un noyau
qui accepte des modèles d'un noyau qui accepte UN format.

CE QU'IL COUVRE : `<worldbody>` imbriqué, `<default>` par classe avec
`childclass`, `<geom>` de type capsule / sphere / box / cylinder (masse par
densité ou `mass=` explicite), `<joint>` hinge et slide, `<freejoint>`.
Plusieurs charnières concourantes sur un corps deviennent la liaison qui
libère leurs axes.

CE QU'IL NE COUVRE PAS, déclaré : `mesh` et `hfield` (géométrie externe),
`tendon`, `equality`, `actuator`, les `<geom>` de type
`plane` (pris comme bâti, sans masse) et les orientations par `xyaxes` /
`zaxis` / `axisangle` — refusées explicitement plutôt que devinées.

Les modèles du **MuJoCo Menagerie** entrent par ce chemin : leur géométrie est
en `mesh` (fichiers externes qu'on ne lit pas) mais ils portent tous leur
`<inertial>`. Sans la priorité donnée à celui-ci, un bras industriel entier
serait sans masse alors que le fichier la porte.
"""

import xml.etree.ElementTree as ET

import numpy as np

from ._vinkulum import Noyau

RHO = 1000.0          # masse volumique par défaut de MuJoCo, kg/m³


def _fini(v, quoi):
    """Refuse NaN et inf — un modèle à `nan` kg s'intègre sans erreur en
    produisant du `nan` partout, et c'est le pire des silences."""
    if not np.all(np.isfinite(np.asarray(v, float))):
        raise ValueError(f"{quoi} : valeur non finie ({v})")
    return v


def _quat(q):
    """quaternion MuJoCo (w, x, y, z) → matrice de rotation."""
    w, x, y, z = q
    nn = (w * w + x * x + y * y + z * z) ** 0.5
    if nn < 1e-12:
        return np.eye(3)
    w, x, y, z = w / nn, x / nn, y / nn, z / nn
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def _euler(v, deg):
    """euler MuJoCo, séquence xyz extrinsèque par défaut."""
    a = np.radians(v) if deg else np.asarray(v, float)
    cx, sx, cy, sy, cz, sz = (np.cos(a[0]), np.sin(a[0]), np.cos(a[1]),
                              np.sin(a[1]), np.cos(a[2]), np.sin(a[2]))
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return rz @ ry @ rx


def _rot_de(e, deg):
    """orientation d'un élément : quat, euler ou rien."""
    if e.get("quat") is not None:
        q = _fini(_f(e.get("quat"), 4), "quaternion")
        if np.linalg.norm(q) < 1e-12:
            raise ValueError(f"quaternion NUL ({e.get('quat')}) — pas une rotation")
        return _quat(q)
    if e.get("euler") is not None:
        return _euler(_f(e.get("euler"), 3), deg)
    if e.get("xyaxes") is not None:
        v = _f(e.get("xyaxes"), 6)
        x = np.asarray(v[:3], float)
        x /= np.linalg.norm(x)
        y = np.asarray(v[3:], float) - np.dot(v[3:], x) * x   # orthogonalisé, comme MuJoCo
        y /= np.linalg.norm(y)
        return np.column_stack([x, y, np.cross(x, y)])
    if e.get("zaxis") is not None:
        z = np.asarray(_f(e.get("zaxis"), 3), float)
        z /= np.linalg.norm(z)
        t = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        x = np.cross(t, z)
        x /= np.linalg.norm(x)
        return np.column_stack([x, np.cross(z, x), z])
    if e.get("axisangle") is not None:
        v = _f(e.get("axisangle"), 4)
        ax = np.asarray(v[:3], float)
        nn = np.linalg.norm(ax)
        if nn < 1e-12:
            return np.eye(3)
        ax = ax / nn
        th = np.radians(v[3]) if deg else v[3]
        k = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        return np.eye(3) + np.sin(th) * k + (1 - np.cos(th)) * (k @ k)
    return np.eye(3)


def _f(s, n=None):
    v = [float(x) for x in str(s).split()]
    return v if n is None else (v + [0.0] * n)[:n]


def _inertie_geom(g, rho=RHO):
    """(masse, CdM local, tenseur 3×3 au CdM, repère local) d'un `<geom>`.

    Les formules sont celles des primitives ; la capsule est un cylindre plus
    deux hémisphères, et c'est la seule qui demande un Steiner (le CdM d'un
    hémisphère est à 3r/8 de son plan).
    """
    typ = g.get("type", "capsule")
    siz = _f(g.get("size", "0"))
    ft = g.get("fromto")
    p1 = np.zeros(3)
    p2 = np.zeros(3)
    if ft is not None:
        v = _f(ft, 6)
        p1, p2 = np.array(v[:3]), np.array(v[3:])
    c = (p1 + p2) / 2
    ax = p2 - p1
    lg = float(np.linalg.norm(ax))
    e1 = ax / lg if lg > 1e-12 else np.array([0.0, 0.0, 1.0])
    t = np.array([0.0, 0.0, 1.0]) if abs(e1[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e2 = np.cross(e1, t)
    e2 /= np.linalg.norm(e2)
    R = np.column_stack([e1, e2, np.cross(e1, e2)])   # e1 = axe de la primitive

    if typ == "sphere":
        r = siz[0]
        m = rho * 4 / 3 * np.pi * r ** 3
        j = 0.4 * m * r * r
        return m, c, np.eye(3) * j, np.eye(3)
    if typ == "ellipsoid":
        a_, b_, d_ = (siz + [0, 0, 0])[:3]
        m = rho * 4 / 3 * np.pi * a_ * b_ * d_
        j = np.diag([m * (b_ * b_ + d_ * d_) / 5, m * (a_ * a_ + d_ * d_) / 5,
                     m * (a_ * a_ + b_ * b_) / 5])
        return m, c, j, np.eye(3)
    if typ == "box":
        a, b, d = (siz + [0, 0, 0])[:3]
        m = rho * 8 * a * b * d
        j = np.diag([m * (b * b + d * d) / 3, m * (a * a + d * d) / 3,
                     m * (a * a + b * b) / 3])
        return m, c, j, np.eye(3)
    if typ in ("capsule", "cylinder"):
        r = siz[0]
        h = lg if ft is not None else 2 * siz[1]
        mc = rho * np.pi * r * r * h
        ja = 0.5 * mc * r * r
        jt = mc * (3 * r * r + h * h) / 12
        if typ == "capsule":
            ms = rho * 2 / 3 * np.pi * r ** 3                # un hémisphère
            ja += 2 * 0.4 * ms * r * r
            # hémisphère : I transverse au CdM propre, puis Steiner vers le centre
            jh = 0.4 * ms * r * r - ms * (3 * r / 8) ** 2
            d = h / 2 + 3 * r / 8
            jt += 2 * (jh + ms * d * d)
            mc += 2 * ms
        return mc, c, np.diag([ja, jt, jt]), R
    raise ValueError(f"geom de type « {typ} » non couvert")


def _defauts(root):
    """Classes de `<default>`, aplaties avec héritage du parent."""
    out = {}

    def rec(e, herite):
        cl = e.get("class", "main")
        cur = {k: dict(v) for k, v in herite.items()}
        for ch in e:
            if ch.tag == "default":
                continue
            cur.setdefault(ch.tag, {}).update(ch.attrib)
        out[cl] = cur
        for ch in e.findall("default"):
            rec(ch, cur)

    d = root.find("default")
    if d is not None:
        rec(d, {})
    return out


def charge(chemin, g=(0.0, 0.0, -9.81), rho=RHO):
    """Construit un `Noyau` depuis un fichier MJCF. Rend (noyau, infos)."""
    try:
        root = ET.parse(chemin).getroot()
    except ET.ParseError as ex:
        raise ValueError(f"{chemin} : XML invalide — {ex}") from ex
    if root.tag != "mujoco":
        raise ValueError(f"{chemin} : racine <{root.tag}>, attendu <mujoco>")
    cp = root.find("compiler")
    deg = (cp is None) or cp.get("angle", "degree") == "degree"
    dfl = _defauts(root)
    wb = root.find("worldbody")
    if wb is None:
        raise ValueError(f"{chemin} : pas de <worldbody>")
    n = Noyau(list(g))
    idx, notes, njoints, mtot, racines_libres = {}, [], 0, [0.0], [0]
    # ⚠ UN <include> NON RÉSOLU peut rendre le modèle INCOMPLET — `aloha.xml`
    # inclut ses actionneurs et sa seconde moitié. Mais l'acrobot inclut des
    # MATÉRIAUX (`common/visual.xml`), et refuser l'aurait rejeté alors qu'il
    # se charge parfaitement. On ne peut pas savoir sans lire le fichier
    # inclus : donc on AVERTIT et on laisse l'appelant juger. Un garde-fou qui
    # crie au loup finit ignoré — c'est écrit dans ce dépôt depuis le premier
    # jour, et mon refus criait au loup sur l'acrobot.
    inc = [e.get("file") for e in root.iter("include")]
    if inc:
        notes.append(f"{len(inc)} <include> non résolu(s) ({', '.join(inc[:2])}) — "
                     "le modèle PEUT être incomplet")

    def attrs(e, cls):
        """attributs de `e` complétés par sa classe de défauts"""
        base = dict(dfl.get(e.get("class", cls), {}).get(e.tag, {}))
        base.update(e.attrib)
        return base

    def rec(b, parent, pos_p, rot_p, cls):
        nonlocal njoints
        cls = b.get("childclass", cls)
        rot = rot_p @ _rot_de(b, deg)
        pos = pos_p + rot_p @ np.array(_f(b.get("pos", "0 0 0"), 3))
        # ── masse et inertie ──
        # <inertial> a PRIORITÉ : c'est ce que donnent les modèles du
        # Menagerie, dont la géométrie est en mesh (fichiers externes qu'on ne
        # lit pas). Sans cette priorité, un bras industriel entier est sans
        # masse alors que le fichier la porte.
        ine = b.find("inertial")
        if ine is not None:
            mt = _fini(float(ine.get("mass")), f"masse de « {b.get('name')} »")
            if mt < 0.0:
                raise ValueError(f"masse NEGATIVE ({mt}) sur « {b.get('name')} »")
            ct = np.array(_f(ine.get("pos", "0 0 0"), 3))
            ri = _rot_de(ine, deg)
            if ine.get("diaginertia") is not None:
                jd = np.diag(_f(ine.get("diaginertia"), 3))
            elif ine.get("fullinertia") is not None:
                v = _f(ine.get("fullinertia"), 6)
                jd = np.array([[v[0], v[3], v[4]], [v[3], v[1], v[5]],
                               [v[4], v[5], v[2]]])
            else:
                raise ValueError(f"<inertial> de « {b.get('name')} » sans inertie")
            _fini(jd.ravel(), f"inertie de « {b.get('name')} »")
            jt = ri @ jd @ ri.T
        else:
            mt, ct, jt = 0.0, np.zeros(3), np.zeros((3, 3))
            for gm in b.findall("geom"):
                a = attrs(gm, cls)
                if a.get("type") in ("plane", "hfield", "sdf"):
                    continue
                if a.get("type") == "mesh":
                    # ⚠ UN MESH PEUT PORTER SA MASSE. On ne lit pas la
                    # géométrie externe, mais `mass=` est là et l'ignorer
                    # jette la masse du corps : le Stretch de Hello Robot
                    # sortait à 0,050 kg au lieu de ~25. On prend la masse et
                    # on DÉCLARE que l'inertie de rotation, elle, est inconnue.
                    mv = float(a.get("mass", 0.0))
                    if mv <= 0.0:
                        continue
                    jmin = np.eye(3) * max(mv * 1e-4, 1e-12)
                    ct = (ct * mt + np.zeros(3) * mv) / (mt + mv) if mt + mv > 0 else ct
                    jt = jt + jmin
                    mt += mv
                    notes.append(f"{b.get('name')} : geom `mesh` de {mv:.4f} kg — masse "
                                 f"prise, inertie de ROTATION inconnue (pas de maillage lu)")
                    continue
                e2 = ET.Element("geom", a)
                m_, c_, j_, R_ = _inertie_geom(e2, rho)
                if a.get("mass") is not None:
                    mv = float(a["mass"])
                    j_ = j_ * (mv / m_ if m_ > 0 else 1.0)
                    m_ = mv
                jw = R_ @ j_ @ R_.T
                ct = (ct * mt + c_ * m_) / (mt + m_) if mt + m_ > 0 else ct
                jt = jt + jw
                mt += m_
        if mt <= 0.0:
            mt = 1e-9
            jt = np.eye(3) * 1e-12
            notes.append(f"{b.get('name')} : aucun geom massique — masse minimale imposée")
        mtot[0] += mt
        i = n.corps(b.get("name") or f"b{len(idx)}", mt, list(jt.flatten()),
                    list(pos + rot @ ct), rot=list(rot.flatten()))
        idx[b.get("name")] = i

        jts = [j for j in b.findall("joint")]
        libre = b.find("freejoint") is not None or any(
            attrs(j, cls).get("type") == "free" for j in jts)
        if libre and parent is None:
            racines_libres[0] += 1
        if not libre:
            typs = [attrs(j, cls).get("type", "hinge") for j in jts]
            if any(t not in ("hinge", "slide", "ball") for t in typs):
                raise ValueError(f"body « {b.get('name')} » : type de joint "
                                 f"{set(typs)} non couvert")
            pj = (pos + rot @ np.array(_f(attrs(jts[0], cls).get("pos", "0 0 0"), 3))
                  if jts else pos)
            axes = [rot @ np.array(_f(attrs(j, cls).get("axis", "0 0 1"), 3)) for j in jts]
            # le repère de liaison porte les axes LIBÉRÉS en premier
            br, bt = [0, 1, 2], [0, 1, 2]
            ra = np.eye(3)
            # `ball` = ROTULE, et le noyau la fait depuis toujours (bloque_r
            # vide). La refuser était un refus injustifié : la brique existe.
            if "ball" in typs:
                br = []
                if any(t == "slide" for t in typs):
                    raise ValueError(f"body « {b.get('name')} » : ball + slide "
                                     "non couvert")
            elif len(jts) == 1:
                a0 = axes[0] / max(np.linalg.norm(axes[0]), 1e-30)
                t_ = np.array([0.0, 0.0, 1.0]) if abs(a0[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
                e2_ = np.cross(a0, t_)
                e2_ /= np.linalg.norm(e2_)
                ra = np.column_stack([a0, e2_, np.cross(a0, e2_)])
                (br if typs[0] == "hinge" else bt).remove(0)
            else:
                for t, ax in zip(typs, axes):
                    k = int(np.argmax(np.abs(ax)))
                    cible = br if t == "hinge" else bt
                    if k in cible:
                        cible.remove(k)
            # `pa` et `ra` sont dans le repère du CORPS parent
            rp = rot_p if parent is not None else np.eye(3)
            org = pos_p if parent is not None else np.zeros(3)
            n.liaison(f"j{len(idx)}", parent, i, pa=list(rp.T @ (pj - org)),
                      ra=list((rp.T @ ra).flatten()), bloque_t=bt, bloque_r=br)
            njoints += 1
        for ch in b.findall("body"):
            rec(ch, i, pos, rot, cls)

    for b in wb.findall("body"):
        rec(b, None, np.zeros(3), np.eye(3), "main")
    return n, dict(noms=idx, notes=notes, n_corps=len(idx), n_liaisons=njoints,
                   masse=mtot[0], libre=racines_libres[0] > 0)


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
    """Autotest sur deux modèles MJCF publics, pris tels quels.

    DEUX CONTRÔLES QUI MORDENT, parce que « ça tourne sans planter » ne prouve
    rien :

    · le humanoid porte un `<freejoint>` — rien ne le retient — et ses masses
      sont calculées depuis ses capsules. Lâché, il doit tomber EXACTEMENT de
      g·t²/2 ; c'est ce qui valide à la fois la géométrie massique et le
      degré de liberté libre ;
    · l'acrobot a deux charnières d'axe y. Un couple autour de y doit le faire
      bouger, un couple autour de x ou de z ne doit RIEN faire. Sans ce
      contrôle, une liaison qui bloquerait tout passerait le banc.
    """
    import os
    d = os.path.join(os.path.dirname(__file__), "donnees")
    print("╔═ vinkulum — MJCF : un second format, de sémantique différente")
    out = {}
    for f, nc, mref in (("acrobot.xml", 2, None), ("humanoid.xml", 16, 40.0),
                        ("ur5e.xml", 7, 20.6), ("go1.xml", 13, 12.0)):
        n, info = charge(os.path.join(d, f))
        phi = max(abs(x) for x in n.phi()) if len(n.phi()) else 0.0
        n.simule(0.1, 1e-3, tous=10 ** 9)
        e = n.etat()[1]
        mtot = sum(1 for _ in info["noms"])
        ref = "" if mref is None else f" · publie ~{mref:.1f}"
        print(f"║ {f:14s} {info['n_corps']:2d} corps · {info['n_liaisons']:2d} liaisons · "
              f"{info['masse']:7.3f} kg{ref} · |Φ| {phi:.0e}")
        for s in info["notes"][:2]:
            print(f"║     {s}")
        assert info["n_corps"] == nc, (f, info["n_corps"])
        assert all(np.all(np.isfinite(p)) for p in e), ("etat non fini", f)
        # · la masse doit tomber sur ce que le constructeur (ou MuJoCo) publie.
        #   Pour le humanoid elle est CALCULEE depuis les capsules : c'est la
        #   formule d'inertie qui est jugee, pas la lecture d'un attribut.
        if mref is not None:
            assert abs(info["masse"] - mref) / mref < 0.10, ("masse hors bande", f,
                                                             info["masse"], mref)
        out[f] = info

    # le humanoid est LIBRE : il doit tomber exactement de g t²/2
    nh, _ = charge(os.path.join(d, "humanoid.xml"))
    p0 = np.array(nh.etat()[1])
    nh.simule(0.2, 1e-4, tous=10 ** 9)
    dz = (np.array(nh.etat()[1]) - p0)[:, 2]
    chute = float(np.max(np.abs(dz - (-9.81 * 0.04 / 2))))
    print(f"║ humanoid lache (freejoint) : chute {dz.mean():+.6f} m en 0,2 s · "
          f"g·t²/2 = {-9.81 * 0.04 / 2:+.6f} · ecart {chute:.1e}")
    assert chute < 1e-6, ("la chute libre ne suit pas g — masse ou freejoint faux", chute)

    # l'acrobot : le couple ne passe QUE par l'axe declare
    rep = {}
    for ax, nom in (([0.0, 1.0, 0.0], "y"), ([1.0, 0.0, 0.0], "x"), ([0.0, 0.0, 1.0], "z")):
        na, _ = charge(os.path.join(d, "acrobot.xml"))
        q0 = np.array(na.etat()[1])
        na.effort(0, [0.0, 0.0, 0.0], list(np.asarray(ax) * 20.0))
        na.simule(0.3, 1e-4, tous=10 ** 9)
        rep[nom] = float(np.max(np.abs(np.array(na.etat()[1]) - q0)))
    print(f"║ acrobot, couple de 20 N·m : autour de y (l'axe du joint) {rep['y']:.5f} m · "
          f"autour de x {rep['x']:.0e} · de z {rep['z']:.0e}")
    assert rep["y"] > 0.1, ("l'axe du joint ne libere rien", rep)
    assert rep["x"] < 1e-9 and rep["z"] < 1e-9, ("un axe cense etre bloque bouge", rep)
    out["controles"] = dict(chute=chute, couple=rep)
    # ── LE LOT, pas quatre cas choisis ──
    # Douze robots du MuJoCo Menagerie, pris tels quels. Ce qui compte ici
    # n'est pas qu'un modèle passe mais le TAUX : choisir ses cas est
    # exactement ce qu'un banc ne doit pas faire, et c'est le reproche qu'on
    # peut adresser à tous les bancs de cette suite pris un par un.
    dm = os.path.join(d, "menagerie")
    lot, echecs, refuses = [], [], []
    for f in sorted(os.listdir(dm)):
        if not f.endswith(".xml"):
            continue
        try:
            nl, il = charge(os.path.join(dm, f))
            nl.simule(0.05, 1e-3, tous=10 ** 9)
            assert all(np.all(np.isfinite(p)) for p in nl.etat()[1])
            # ⚠ LE CONTRÔLE DE CONSERVATION NE VAUT QUE POUR UN SYSTÈME LIBRE.
            # La moitié du Menagerie est BOULONNÉE au sol (les bras
            # industriels) : le bâti y exerce un couple, et L n'a aucune
            # raison d'être conservé — mesuré, 9e-2 sur le Panda contre 1e-10
            # sur les quadrupèdes. Appliquer le critère aux deux aurait fait
            # accuser le solveur d'un fait de modèle.
            dl = None
            if il["libre"]:
                nl2, _ = charge(os.path.join(dm, f), g=(0, 0, 0))
                dl = _conservation(nl2)
            lot.append((f, il["n_corps"], il["masse"], dl))
        except ValueError as ex:                     # refus MOTIVÉ du lecteur
            refuses.append((f, str(ex)[:66]))
        except Exception as ex:                      # noqa: BLE001 — on VEUT tout attraper
            echecs.append((f, f"{type(ex).__name__}: {ex}"))
    print(f"║ ── le LOT : {len(lot)} robots du MuJoCo Menagerie charges, "
          f"{len(refuses)} refuses avec motif ──")
    for k in range(0, len(lot), 6):
        print("║   " + " · ".join(f"{f[:-4]} {c}c {m:.1f}kg" for f, c, m, _ in lot[k:k + 6]))
    for f, e in refuses:
        print(f"║   refus {f} : {e}")
    lib = [x for x in lot if x[3] is not None]
    print(f"║   conservation du moment cinetique : dL/L max {max(x[3] for x in lib):.1e} "
          f"sur les {len(lib)} modeles LIBRES ({len(lot) - len(lib)} sont boulonnes au sol,")
    print("║   le critere ne s'y applique pas — le bati y exerce un couple)")
    assert max(x[3] for x in lib) < 1e-6, ("un modele libre ne conserve pas L", lib)
    assert len(lib) >= 5, ("plus assez de modeles libres pour juger", lib)
    for f, e in echecs:
        print(f"║   ECHEC {f} : {e[:70]}")
    assert not echecs, ("des modeles publics ne passent plus", echecs)
    assert len(lot) >= 19, ("le lot a maigri", len(lot))
    out["menagerie"] = lot

    print("╚ MJCF OK — deux formats de sémantiques différentes lus par le même noyau")
    return out


if __name__ == "__main__":
    demo()

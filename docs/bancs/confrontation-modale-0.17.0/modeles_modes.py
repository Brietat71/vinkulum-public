"""Modèles modaux explicites et références mécaniques indépendantes."""
import numpy as np


def reference_console(nb,combien=6):
    """Flexibilités exactes du modèle linéaire discrétisé, assemblage indépendant.

    v_i=l Σγ_j+l² Σ(i-j+1/2)κ_j ; θ_i=l Σκ_j.
    Les grandes valeurs propres de sqrt(M) K⁻¹ sqrt(M) évitent une
    soustraction au bas du spectre de K. Aucun appel au noyau.
    """
    from scipy.linalg import eigvalsh
    E,rho,L,b,h=70e9,2700.,1.,.02,.005
    area=b*h;Iy=b*h**3/12;Iz=h*b**3/12;J=Iy+Iz;dl=L/nb
    ga=E/2.6*area;ei=E*Iy
    i,j=np.meshgrid(np.arange(1,nb+1,dtype=np.int64),np.arange(1,nb+1,dtype=np.int64),indexing='ij')
    a=np.minimum(i,j)
    rr4=a*(2*i+1)*(2*j+1)-2*(i+j+1)*a*(a+1)+2*a*(a+1)*(2*a+1)//3
    rl2=a*(2*i-a)
    vv=dl/ga*a+dl**3/ei*(rr4/4.)
    vt=dl**2/ei*(rl2/2.)
    tt=dl/ei*a
    compliance=np.block([[vv,vt],[vt.T,tt]])
    weight=np.ones(nb);weight[-1]=.5
    values=[]
    for moment in (Iy,Iz):
        root=np.sqrt(np.concatenate((rho*area*dl*weight,rho*moment*dl*weight)))
        inverse=root[:,None]*compliance*root[None,:]
        eig=eigvalsh(inverse,subset_by_index=(2*nb-min(combien,2*nb),2*nb-1))
        values.extend((1./eig).tolist())
    k=np.arange(1,min(combien,nb)+1)
    for wave in (E/rho,E/2.6/rho):
        values.extend((4*wave/dl**2*np.sin((2*k-1)*np.pi/(4*nb))**2).tolist())
    return np.sort(values)[:combien]


def reference_chaine(nb,gravite=False):
    """Angles absolus des barres ; énergie cinétique et potentiel analytiques."""
    from scipy.linalg import eigvalsh
    length=mass=1./nb
    i,j=np.indices((nb,nb))
    # 6 M/(m l²) a des coefficients entiers : conserver ces sommes
    # exactement avant l'unique remise à l'échelle en flottants.
    integers=6*(nb-np.maximum(i,j))-3-np.eye(nb,dtype=np.int64)
    if gravite:
        m=integers*(mass*length**2/6.)
        root=np.sqrt(9.81*mass*length*(nb-np.arange(nb)-.5))
        inverse=m/root[:,None]/root[None,:]
    else:
        # Angles relatifs : K = k I et M = Tᵀ M_abs T. Les PLUS
        # GRANDES valeurs de M/k donnent les petites pulsations par inverse,
        # sans les chercher au bas d'un spectre mal conditionné.
        relative=integers[::-1,::-1].cumsum(axis=0).cumsum(axis=1)[::-1,::-1]
        inverse=relative*(mass*length**2/(60.*nb))
    return 1./eigvalsh(inverse)[::-1]


def vinkulum(family,nb):
    from vinkulum import Noyau
    if family=='redondante':
        from vinkulum.test_noyau import _cascade_initiale
        return _cascade_initiale(nb,np.eye(3))
    if family=='console':
        E,rho,L,b,h=70e9,2700.,1.,.02,.005
        area=b*h;Iy=b*h**3/12;Iz=h*b**3/12;J=Iy+Iz;dl=L/nb
        n=Noyau([0,0,0])
        for i in range(nb+1):
            weight=.5 if i in (0,nb) else 1.
            mass=rho*area*dl*weight
            inertia=rho*dl*weight*np.diag([J,Iy,Iz])
            n.corps(str(i),mass,inertia.ravel().tolist(),[i*dl,0,0])
        n.liaison('encastrement',None,0)
        for i in range(nb):n.poutre(str(i),i,i+1,E*area,E/2.6*area,E/2.6*J,E*Iy)
        return n
    if family not in ('gravite','articulee'):raise ValueError('famille inconnue')
    length=mass=1./nb
    n=Noyau([0,0,-9.81 if family=='gravite' else 0.])
    for i in range(nb):
        n.corps(str(i),mass,(np.eye(3)*mass*length**2/12).ravel().tolist(),[0,0,-(i+.5)*length])
        a=None if not i else i-1
        n.liaison(str(i),a,i,pa=[0,0,0 if not i else -length/2],bloque_r=[0,2])
        if family=='articulee':n.couple('rappel'+str(i),a,i,[0,1,0],('ressort',[10.*nb,0.,0.]))
    return n


def exudyn(nb,arbre=False):
    """Même chaîne à ressorts, construite par l'interface officielle Exudyn."""
    import exudyn as exu
    from exudyn.itemInterface import NodeRigidBodyEP,ObjectRigidBody,ObjectGround,MarkerBodyRigid,ObjectJointGeneric,ObjectConnectorRigidBodySpringDamper
    sc=exu.SystemContainer();s=sc.AddSystem();length=mass=1./nb
    if arbre:
        from exudyn.itemInterface import ObjectKinematicTree,NodeGenericODE2
        node=s.AddNode(NodeGenericODE2(referenceCoordinates=[0.]*nb,initialCoordinates=[0.]*nb,
                                     initialCoordinates_t=[0.]*nb,numberOfODE2Coordinates=nb))
        s.AddObject(ObjectKinematicTree(nodeNumber=node,jointTypes=[exu.JointType.RevoluteY]*nb,
            linkParents=list(range(-1,nb-1)),jointTransformations=exu.Matrix3DList([np.eye(3)]*nb),
            jointOffsets=exu.Vector3DList([[0,0,0]]+[[0,0,-length]]*(nb-1)),
            linkInertiasCOM=exu.Matrix3DList([np.eye(3)*mass*length**2/12]*nb),
            linkCOMs=exu.Vector3DList([[0,0,-length/2]]*nb),linkMasses=[mass]*nb,
            jointPositionOffsetVector=[0.]*nb,jointPControlVector=[10.*nb]*nb))
        s.Assemble()
        return sc,s
    ground=s.AddObject(ObjectGround())
    bodies=[]
    for i in range(nb):
        node=s.AddNode(NodeRigidBodyEP(referenceCoordinates=[0,0,-(i+.5)*length,1,0,0,0]))
        body=s.AddObject(ObjectRigidBody(nodeNumber=node,physicsMass=mass,
                       physicsInertia=[mass*length**2/12]*3+[0,0,0]))
        first=s.AddMarker(MarkerBodyRigid(bodyNumber=ground if not i else bodies[-1],
                         localPosition=[0,0,0 if not i else -length/2]))
        second=s.AddMarker(MarkerBodyRigid(bodyNumber=body,localPosition=[0,0,length/2]))
        s.AddObject(ObjectJointGeneric(markerNumbers=[first,second],constrainedAxes=[1,1,1,1,0,1]))
        stiff=np.zeros((6,6));stiff[4,4]=10.*nb
        s.AddObject(ObjectConnectorRigidBodySpringDamper(markerNumbers=[first,second],stiffness=stiff))
        bodies.append(body)
    s.Assemble()
    return sc,s

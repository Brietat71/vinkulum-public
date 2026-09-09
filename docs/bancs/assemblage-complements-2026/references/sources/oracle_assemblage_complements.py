"""Oracle Decimal : Dirichlet par branches, puis équilibre des ports communs.

Aucune base réduite ni Schur candidat n'entre dans cette référence.
Les pivots intérieurs par blocs doivent être inversibles ; le Schur d'une
pièce isolée peut être singulier. Convergence de précision, pas certification.
"""
from decimal import localcontext
import numpy as np
from oracle_champs_couples import OracleChampsCouples
from reference_ports_precision import _decimal,_ZERO,_UN,_resout_bloc,_resultat,_precision


def produit(a,b):
    return [[sum((x*y for x,y in zip(row,col,strict=True)),_ZERO) for col in zip(*b)] for row in a]

def transpose(a):return list(map(list,zip(*a)))

def decimal_matrix(a):
    if np.iscomplexobj(a):raise ValueError('matrice réelle requise')
    a=np.asarray(a)
    if a.ndim!=2:raise ValueError('matrice requise')
    return [[_decimal(v,'coefficient') for v in row] for row in a]


class OracleDirichlet(OracleChampsCouples):
    def _preparer_frequence(self,omega):
        if np.iscomplexobj(omega):raise ValueError('fréquence réelle requise')
        omega=_decimal(omega,'omega')
        if omega<0:raise ValueError('fréquence non négative requise')
        if omega==self._omega:return
        z,p=omega*omega,self.p
        def diagonal(i):
            return [[self._diagonaux[i][j][k]-z*self._masse_diagonaux[i][j][k] for k in range(p)] for j in range(p)]
        pivot=diagonal(0);transfers=[]
        for i,k in enumerate(self._liaisons):
            link=[[k[j][l]-z*self._masse_liaisons[i][j][l] for l in range(p)] for j in range(p)]
            solved=_resout_bloc(pivot,link);next_block=diagonal(i+1)
            for j in range(p):
                for l in range(j,p):
                    jl=sum((link[h][j]*solved[h][l] for h in range(p)),_ZERO)
                    lj=sum((link[h][l]*solved[h][j] for h in range(p)),_ZERO)
                    next_block[j][l]-=(jl+lj)/2;next_block[l][j]=next_block[j][l]
            transfers.append(solved);pivot=next_block
        # Le Schur terminal n'est pas inversé : il peut être singulier.
        self._schur,self._transferts,self._omega=pivot,transfers,omega
        self.compteurs['condensations_frequence']+=1
        self.compteurs['resolutions_blocs']+=len(transfers)
    def schur(self,omega,retour_decimal=False):
        with localcontext(self._contexte):
            self._preparer_frequence(omega);return _resultat(self._schur,retour_decimal)
    def champ_dirichlet(self,omega,port,retour_decimal=False):
        with localcontext(self._contexte):
            self._preparer_frequence(omega);current=decimal_matrix(port)
            if len(current)!=self.p or not current[0]:raise ValueError('déplacements de port incompatibles')
            blocks=[current]
            for transfer in reversed(self._transferts):
                current=[[-v for v in row] for row in produit(transfer,current)];blocks.append(current)
            return _resultat([row for block in reversed(blocks) for row in block],retour_decimal)
    def reponse(self,*args,**kwargs):raise ValueError('utiliser champ_dirichlet avec des déplacements prescrits')


class OracleAssemblage:
    def __init__(self,pieces,applications,facteur_externe,masse_externe,*,dps=70):
        if not pieces or len(pieces)!=len(applications):raise ValueError('pièces et applications appariées requises')
        self.context=_precision(dps);self.g=np.asarray(applications[0]).shape[1]
        self.apps=[decimal_matrix(a) for a in applications]
        self.pieces=[OracleDirichlet(d,m,p=len(a),dps=dps) for (d,m),a in zip(pieces,self.apps,strict=True)]
        if any(len(a[0])!=self.g for a in self.apps):raise ValueError('dimension globale incohérente')
        de=decimal_matrix(facteur_externe);self.me=decimal_matrix(masse_externe)
        if len(self.me)!=self.g or any(len(row)!=self.g for row in self.me):raise ValueError('masse externe incompatible')
        with localcontext(self.context):self.ke=produit(transpose(de),de)
    def reponses(self,omega,forces,retour_decimal=False):
        with localcontext(self.context):
            if np.iscomplexobj(omega):raise ValueError('fréquence réelle requise')
            w=_decimal(omega,'omega');f=decimal_matrix(forces)
            if w<0 or len(f)!=self.g or not f[0]:raise ValueError('fréquence ou forces incompatibles')
            s=[[self.ke[i][j]-w*w*self.me[i][j] for j in range(self.g)] for i in range(self.g)]
            for piece,a in zip(self.pieces,self.apps,strict=True):
                local=piece.schur(w,retour_decimal=True);added=produit(transpose(a),produit(local,a))
                for i in range(self.g):
                    for j in range(self.g):s[i][j]+=added[i][j]
            identity=[[_UN if i==j else _ZERO for j in range(self.g)] for i in range(self.g)]
            u=produit(_resout_bloc(s,identity),f)
            fields=[piece.champ_dirichlet(w,produit(a,u),retour_decimal=retour_decimal) for piece,a in zip(self.pieces,self.apps,strict=True)]
            return dict(champs=fields,champ_ports=_resultat(u,retour_decimal))

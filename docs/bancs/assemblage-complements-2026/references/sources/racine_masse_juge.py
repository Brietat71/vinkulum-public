"""Facteur massique indépendant du candidat, avec défaut exact Fraction.

La bande limite le stockage du juge, pas les modèles acceptés par le
contrôleur. Le facteur flottant n'est jamais supposé exact : son défaut
est majoré relativement à diag(M), puis à M via une comparaison prouvée.
"""
from decimal import Decimal
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import cholesky_banded
from scipy.sparse import csr_matrix

from comparaison_masse import certifier_comparaison_masse
from matrices_certificat import matrice_entree
from inertie_complement_dirigee import IntervallesSignes
from vinkulum._ports.certificat_spectral import _float_dirige


class RacineMasseJuge:
    def __init__(self,m,*,bibliotheque,budget_coefficients=2_000_000,budget_produits=2_000_000):
        m=matrice_entree(m,'M');n=m.shape[0]
        if m.shape!=(n,n) or n==0 or np.any((m-m.T).data!=0):raise ValueError('masse carrée exactement symétrique requise')
        self.comparaison=certifier_comparaison_masse(m,.25,2.,bibliotheque=bibliotheque)
        coo=m.tocoo();width=int(np.max(np.abs(coo.row-coo.col),initial=0))
        if n*(width+1)>budget_coefficients:raise ValueError('budget de stockage bande du juge dépassé')
        band=np.zeros((width+1,n));keep=coo.row>=coo.col
        band[coo.row[keep]-coo.col[keep],coo.col[keep]]=coo.data[keep]
        lower=cholesky_banded(band,lower=True,check_finite=False)
        rows=[];cols=[];vals=[]
        for k in range(width+1):
            for j in range(n-k):
                v=lower[k,j]
                if v!=0:rows.append(j);cols.append(j+k);vals.append(v)
        if not np.all(np.isfinite(vals)):raise ArithmeticError('facteur de masse non fini')
        self.r=csr_matrix((vals,(rows,cols)),shape=(n,n))
        defect={}
        for i,j,v in zip(coo.row,coo.col,coo.data,strict=True):
            if i<=j:defect[(int(i),int(j))]=-F(float(v))
        products=0
        for row in range(n):
            a,b=self.r.indptr[row:row+2];indices=self.r.indices[a:b];values=[F(float(x)) for x in self.r.data[a:b]]
            for ii,i in enumerate(indices):
                for jj in range(ii,len(indices)):
                    products+=1
                    if products>budget_produits:raise ValueError('budget de produits rationnels du juge dépassé')
                    key=(int(i),int(indices[jj]));defect[key]=defect.get(key,F(0))+values[ii]*values[jj]
        diagonal=list(map(lambda v:F(float(v)),m.diagonal()))
        exponents=[-(math.frexp(float(x))[1]//2) for x in m.diagonal()]
        scales=[F(2)**e for e in exponents]
        sums=[F(0)]*n
        for (i,j),v in defect.items():
            value=abs(v)*scales[i]*scales[j];sums[i]+=value
            if i!=j:sums[j]+=value
        minimum=min(v*s*s for v,s in zip(diagonal,scales,strict=True))
        rho=max(sums)/minimum;alpha=F(self.comparaison['alpha'])
        if rho>=alpha:raise ArithmeticError('défaut du facteur massique hors marge prouvée')
        ratio=(alpha+rho)/(alpha-rho)
        ar=IntervallesSignes(40,1000)
        upper=ar.haut.divide(Decimal(ratio.numerator),Decimal(ratio.denominator))
        self.correction=_float_dirige(ar.racine_superieure(upper),True)
        self.audit=dict(produits_rationnels=products,largeur_bande=width,nnz_facteur=self.r.nnz,
            rho_numerateur=str(rho.numerator),rho_denominateur=str(rho.denominator),
            alpha=self.comparaison['alpha'],correction_relative=self.correction,
            portee='défaut du facteur exact encadré ; calcul des champs pondérés et QR du juge non certifiés')

"""Oracle Decimal de blocs voisins avec masse pleinement couplée.

Le Gram D.T D est assemblé depuis D exact. La masse complète stockée est
soustraite à chaque fréquence, diagonaux de blocs et liaisons compris.
Les références 70/90 chiffres sont convergées, pas certifiées par intervalles.
"""
from decimal import localcontext
import time
import numpy as np

from oracle_champs_ports import OracleChamps,_csr_binary64
from reference_ports_precision import _decimal,_ZERO,_UN,_resout_bloc,_zeros


class OracleChampsCouples(OracleChamps):
    """Même élimination par blocs <=6, sans restriction de masse diagonale.

    D et M peuvent relier deux blocs voisins. Aucun couplage non nul n'est
    ignoré ; une liaison hors bande est refusée. M doit être symétrique mais
    sa positivité ne relève pas de cet oracle d'inversion. Les pivots de
    blocs doivent être inversibles, même si le système global est régulier.
    """
    def __init__(self,d,m,p=6,dps=70):
        start=time.perf_counter()
        m=_csr_binary64(m,'M')
        # Le constructeur parent fournit uniquement l'assemblage du Gram et
        # les compteurs ; sa masse diagonale provisoire n'est jamais utilisée.
        super().__init__(d,np.ones(d.shape[1],dtype=np.float64),p=p,dps=dps)
        if m.shape!=(self.total,self.total) or np.any((m-m.T).data!=0):
            raise ValueError('M exactement symétrique de dimension physique requise')
        self._masse_diagonaux=[_zeros(p) for _ in range(self.nombre_blocs)]
        self._masse_liaisons=[_zeros(p) for _ in range(self.nombre_blocs-1)]
        for i in range(self.total):
            for k in range(m.indptr[i],m.indptr[i+1]):
                j=int(m.indices[k]);v=m.data[k]
                if v==0 or j<i:continue
                bi,ci=divmod(i,p);bj,cj=divmod(j,p)
                if bj-bi>1:raise ValueError('M relie des blocs non voisins')
                target=self._masse_diagonaux[bi] if bi==bj else self._masse_liaisons[bi]
                target[ci][cj]=_decimal(v,'M')
                if bi==bj:target[cj][ci]=target[ci][cj]
        self._masses=None
        self.preparation_s=time.perf_counter()-start

    def _preparer_frequence(self, omega):
        if np.iscomplexobj(omega):
            raise ValueError("fréquence réelle requise")
        omega = _decimal(omega, "omega")
        if omega < 0:
            raise ValueError("omega doit être positif ou nul")
        if omega == self._omega:
            return
        z, p = omega*omega, self.p

        def diagonal(i):
            bloc = [row.copy() for row in self._diagonaux[i]]
            for j in range(p):
                for k in range(p):
                    bloc[j][k] -= z*self._masse_diagonaux[i][j][k]
            return bloc

        pivot = diagonal(0)
        transferts = []
        for i, raideur in enumerate(self._liaisons):
            liaison = [[raideur[j][k]-z*self._masse_liaisons[i][j][k]
                        for k in range(p)] for j in range(p)]
            try:
                resolu = _resout_bloc(pivot, liaison)
            except ValueError as exc:
                raise ValueError(f"pivot intérieur de bloc {i} nul à omega={omega}") from exc
            suivant = diagonal(i+1)
            for j in range(p):
                for k in range(j, p):
                    jk = sum((liaison[l][j]*resolu[l][k] for l in range(p)), _ZERO)
                    kj = sum((liaison[l][k]*resolu[l][j] for l in range(p)), _ZERO)
                    suivant[j][k] -= (jk+kj)/2
                    suivant[k][j] = suivant[j][k]
            transferts.append(resolu)
            pivot = suivant
        identite = [[_UN if i == j else _ZERO for j in range(p)] for i in range(p)]
        try:
            inverse = _resout_bloc(pivot, identite)
        except ValueError as exc:
            raise ValueError(f"Schur terminal singulier à omega={omega}") from exc
        self._transferts, self._schur, self._inverse_schur = transferts, pivot, inverse
        self._omega = omega
        self.compteurs["condensations_frequence"] += 1
        self.compteurs["resolutions_blocs"] += self.nombre_blocs


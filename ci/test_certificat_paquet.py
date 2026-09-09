"""Le pont livré reproduit les preuves archivées et les calculs du prototype."""
import hashlib
from pathlib import Path
import unittest
import numpy as np
from scipy.sparse import csr_matrix
from archive_inertie_binaire import lire_json,preuve_stable
from confronte_ports_exudyn import lire
from vinkulum import _vinkulum
from vinkulum.reduction_contrainte import ReductionContrainte
from vinkulum._ports.inertie_binaire_compile import BibliothequeInertie,certifier_inertie_compilee
from vinkulum._ports.comparaison_masse import BibliothequeComparaison,certifier_comparaison_masse
from condensation_energie import CondensationEnergie
from krylov_contraint import KrylovContraint
from controle_masse_comparee import ControleMasseComparee
from comparaison_masse import BibliothequeComparaison as AncienneBibliotheque,compiler_comparaison
from test_krylov_contraint import donnees
import tempfile

ROOT=Path(__file__).resolve().parent.parent

class PreuvesLivrees(unittest.TestCase):
    def test_sources_cpp_identiques_aux_sources_mesurees(self):
        _,_,inertia,mass=_vinkulum._construction_certificat()
        for text,name in ((inertia,'inertie_binaire_native.cpp'),(mass,'comparaison_masse_native.cpp')):
            self.assertEqual(hashlib.sha256(text.encode()).hexdigest(),hashlib.sha256((ROOT/'ci'/name).read_bytes()).hexdigest())

    def test_trois_preuves_archivees_rejouees_par_extension(self):
        archive=ROOT/'docs/bancs/masse-comparee-2026'
        rows=lire_json(archive,'campagne/bilan.json')
        for n in (32,128,512):
            row=next(r for r in rows if r['n']==n and r['variante']=='masse_comparee')
            d,m,*_=lire(archive/f'references/n{n}-f40.npz');mi=m[:-6,:-6].tocsr()
            with np.load(archive/('campagne/'+row['dossier']+'/bases.npz'),allow_pickle=False) as z:b=z['B']
            expected=row['certificat_complement']
            actual=certifier_inertie_compilee(d[:,:-6],mi,b,expected['lambda_min'],
                permutation=expected['permutation_physique'],bibliotheque=BibliothequeInertie())
            self.assertEqual(preuve_stable(actual),preuve_stable(expected))
            comp=certifier_comparaison_masse(mi,.49,1.51,bibliotheque=BibliothequeComparaison())
            self.assertEqual(preuve_stable(comp),preuve_stable(row['comparaison_masse']))

    def test_algorithmes_historiques_et_publics_identiques(self):
        d,m,i,s,_,phi,_=donnees(8);h=np.diag([4.,.25]);b=m[np.ix_(i,i)]@phi
        public=ReductionContrainte(d,m,i,s,h,.5,gamma=4.,directions_retenues=phi,blocs=1)
        with tempfile.TemporaryDirectory(prefix='vinkulum-equivalence-paquet-') as tmp:
            bib=AncienneBibliotheque(compiler_comparaison(Path(tmp)/'build'))
            old=ControleMasseComparee(KrylovContraint(CondensationEnergie(csr_matrix(d),i,s,h),m,b,phi,4.,.5,blocs=1),bibliotheque=bib)
            for w in (0.,.25,.5):
                a=public.reponses(w,np.eye(2));b=old.reponses(w,np.eye(2))
                for name in ('champ','coordonnees','marge','borne_schur','borne_coordonnees'):
                    np.testing.assert_array_equal(a[name],b[name])
                for name in ('masse','deformation'):
                    for key in ('normes','absolues','relatives'):
                        np.testing.assert_array_equal(a['bornes'][name][key],b['bornes'][name][key])

if __name__=='__main__':unittest.main()

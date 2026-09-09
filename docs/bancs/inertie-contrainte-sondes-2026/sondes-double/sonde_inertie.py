"""Sondes double, sans certificat et sans comparaison chronometrique."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from scipy.linalg import eig_banded, ldl
from scipy.sparse import bmat, coo_matrix, csr_matrix, diags
from scipy.sparse.csgraph import connected_components, reverse_cuthill_mckee

ARCHIVE = Path('/tmp/vinkulum-energie-native-travail/docs/bancs/krylov-contraint-2026/essais')
INPUTS = Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')
OUTPUT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_input(path):
    with np.load(path, allow_pickle=False) as z:
        return [csr_matrix((z[k+'_data'],z[k+'_indices'],z[k+'_indptr']),
                           shape=tuple(z[k+'_shape'])) for k in ('d','m')]

def pattern_stats(a, order):
    c = a[order][:,order].tocoo()
    graph = [set() for _ in order]
    for i,j in zip(c.row,c.col):
        if i != j:
            graph[i].add(int(j))
    degrees = []
    added = 0
    for i in range(len(order)):
        neighbors = sorted(j for j in graph[i] if j > i)
        degrees.append(len(neighbors))
        for pos,j in enumerate(neighbors):
            for k in neighbors[pos+1:]:
                if k not in graph[j]:
                    graph[j].add(k)
                    graph[k].add(j)
                    added += 1
    return dict(nnz=a.nnz, bandwidth=int(np.max(np.abs(c.row-c.col))),
                l_nnz=int(len(order)+sum(degrees)),
                max_future_neighbors=max(degrees), fill_edges=added,
                update_pairs=int(sum(d*(d+1)//2 for d in degrees)))

def scalar_ldl(a):
    """No pivoting; every sign is merely a binary64 observation."""
    c=a.tocoo()
    rows=[{} for _ in range(a.shape[0])]
    for i,j,v in zip(c.row,c.col,c.data):
        if i<=j:
            rows[i][int(j)]=float(v)
    pivots=[]; ratios=[]; rr=[]; cc=[]; vv=[]
    initialmax=float(np.max(np.abs(a.data)))
    workmax=initialmax
    for k,row in enumerate(rows):
        pivot=row.get(k,0.)
        if pivot==0 or not np.isfinite(pivot):
            return dict(refusal='zero or nonfinite pivot', pivot_index=k)
        others=sorted((i,v) for i,v in row.items() if i>k and v!=0.)
        pivots.append(pivot)
        ratios.append(abs(pivot)/max([abs(pivot)]+[abs(v) for _,v in others]))
        rr.append(k);cc.append(k);vv.append(1.)
        for i,v in others:
            rr.append(i);cc.append(k);vv.append(v/pivot)
        for pos,(i,vi) in enumerate(others):
            quotient=vi/pivot
            ri=rows[i]
            for j,vj in others[pos:]:
                result=ri.get(j,0.)-quotient*vj
                ri[j]=result
                workmax=max(workmax,abs(result))
        row.clear()
    p=np.asarray(pivots)
    l=coo_matrix((vv,(rr,cc)),shape=a.shape).tocsr()
    reconstructed=l@diags(p)@l.T
    residual=(reconstructed-a).tocsr()
    residual.eliminate_zeros()
    norminf=lambda x:float(np.max(np.asarray(abs(x).sum(axis=1)).ravel()))
    return dict(inertia=[int(np.count_nonzero(p>0)),int(np.count_nonzero(p<0)),int(np.count_nonzero(p==0))],
                negative_pivot_indices=np.flatnonzero(p<0).tolist(),
                min_abs_pivot=float(np.min(abs(p))), max_abs_pivot=float(np.max(abs(p))),
                final_pivot=float(p[-1]), min_relative_row_pivot=float(min(ratios)),
                largest_l_coefficient=float(max(abs(l.data))),
                active_entry_growth=workmax/initialmax,
                relative_reconstruction_inf=norminf(residual)/norminf(a),
                l_nnz=l.nnz)

def dense_pivoted_ldl(a):
    lu,d,perm=ldl(a.toarray(),lower=True,check_finite=False)
    inertia=[0,0,0]; blocks=[]; i=0
    while i<len(d):
        size=2 if i+1<len(d) and d[i,i+1]!=0 else 1
        ev=np.linalg.eigvalsh(d[i:i+size,i:i+size])
        inertia[0]+=int(np.count_nonzero(ev>0)); inertia[1]+=int(np.count_nonzero(ev<0))
        inertia[2]+=int(np.count_nonzero(ev==0))
        if size==2:
            blocks.append(dict(index=i, original_indices=perm[i:i+2].tolist(), values=d[i:i+2,i:i+2].tolist()))
        i+=size
    return dict(inertia=inertia,two_by_two=blocks,
                moved_coordinates=int(np.count_nonzero(perm!=np.arange(len(perm)))),
                l_numerical_nnz=int(np.count_nonzero(lu)))

def banded_spectrum(a):
    c=a.tocoo(); bw=int(np.max(abs(c.row-c.col)))
    ab=np.zeros((bw+1,a.shape[0]))
    for i,j,v in zip(c.row,c.col,c.data):
        if i>=j:
            ab[i-j,j]=v
    ev=eig_banded(ab,lower=True,eigvals_only=True,check_finite=False)
    return dict(inertia=[int(np.count_nonzero(ev>0)),int(np.count_nonzero(ev<0)),int(np.count_nonzero(ev==0))],
                min_abs_eigenvalue=float(np.min(abs(ev))),
                smallest_eigenvalues=ev[:5].tolist())

def case(n):
    started=time.perf_counter()
    entry=INPUTS/f'n{n}-f40.npz'
    basis=ARCHIVE/f'n{n}-f40-krylov_controle-passage0/bases.npz'
    d,m=read_input(entry)
    with np.load(basis,allow_pickle=False) as z:
        ii=z['indices_interieur']; b=z['B'].copy()
    di=d[:,ii]; mi=m[ii][:,ii].tocsr(); k=(di.T@di).tocsr()
    k.eliminate_zeros(); ni=len(ii)
    structural_di=di.tocsr(copy=True)
    row_counts=np.diff(structural_di.indptr)
    structural_di.data=np.ones(structural_di.nnz)
    structural_k=(structural_di.T@structural_di).tocsr()
    structural_k.data=np.ones(structural_k.nnz)
    structural_c=bmat([[structural_k,csr_matrix(b)],[csr_matrix(b.T),None]],format='csr')
    nc,labels=connected_components(k)
    physical_orders={'natural':np.arange(ni),
                     'physical_rcm':reverse_cuthill_mckee(k,symmetric_mode=True)}
    border=np.array([ni])
    kt=bmat([[k,csr_matrix(b)],[csr_matrix(b.T),None]],format='csr')
    structure={name:pattern_stats(kt,np.concatenate((p,border))) for name,p in physical_orders.items()}
    for name,p in physical_orders.items():
        structure[name]['physical']=pattern_stats(k,p)
    # Global RCM sees the dense multiplier. Stop symbolic expansion if wide.
    global_order=reverse_cuthill_mckee(kt,symmetric_mode=True)
    multiplier_position=int(np.flatnonzero(global_order==ni)[0])
    remaining=global_order[multiplier_position+1:]
    connected_remaining=int(np.count_nonzero(b[remaining[remaining<ni]]))
    structure['global_rcm']=dict(**pattern_stats(kt,global_order),multiplier_position=multiplier_position,
                                multiplier_future_neighbors=connected_remaining,
                                clique_after_multiplier_lower_bound=connected_remaining*(connected_remaining-1)//2)
    fresh=[]
    for path in sorted(ARCHIVE.glob(f'n{n}-f40-krylov_*-passage*/bases.npz')):
        with np.load(path,allow_pickle=False) as z:
            other=z['B']
        fresh.append(dict(path=str(path),sha256=sha(path),nonzero=int(np.count_nonzero(other)),
                          absolute_correlation=float(abs(np.vdot(other,b))/(np.linalg.norm(other)*np.linalg.norm(b)))))
    result=dict(n=n,ni=ni,entry=str(entry),entry_sha256=sha(entry),basis=str(basis),basis_sha256=sha(basis),
                d_nnz=di.nnz,k_nnz=k.nnz,m_nnz=mi.nnz,b_nonzero=int(np.count_nonzero(b)),
                b_abs_min=float(np.min(abs(b))),b_abs_max=float(np.max(abs(b))),
                k_diagonal_range=[float(k.diagonal().min()),float(k.diagonal().max())],
                physical_components=np.bincount(labels).tolist(),structure=structure,fresh_bases=fresh,points=[],
                structural_no_cancellation=dict(k_nnz=structural_k.nnz,
                    d_max_row_nnz=int(max(row_counts)),
                    d_upper_assembly_products=int(np.sum(row_counts*(row_counts+1)//2)),
                    natural=pattern_stats(structural_c,np.arange(ni+1)),
                    physical_rcm=pattern_stats(structural_c,np.r_[reverse_cuthill_mckee(structural_k,symmetric_mode=True),ni])))
    physical_scale=1/np.sqrt(k.diagonal())
    border_scale=1/np.linalg.norm(physical_scale[:,None]*b)
    scaling=np.concatenate((physical_scale,[border_scale]))
    for hz in [60,80,50,40,70,100]:
        gamma=float((2*np.pi*hz)**2)
        a=k-gamma*mi
        c=bmat([[a,csr_matrix(b)],[csr_matrix(b.T),None]],format='csr')
        scaled=(diags(scaling)@c@diags(scaling)).tocsr()
        point=dict(hz=hz,gamma_binary64=gamma,border_scale=float(border_scale),orders={})
        for name,p in physical_orders.items():
            order=np.concatenate((p,border)); cp=scaled[order][:,order]
            probe=scalar_ldl(cp)
            if 'negative_pivot_indices' in probe:
                probe['negative_original_indices']=order[probe['negative_pivot_indices']].tolist()
            if n<=128:
                probe['pivoted_dense_check']=dense_pivoted_ldl(cp)
            point['orders'][name]=probe
        p=physical_orders['physical_rcm']
        physical_scaled=scaled[:ni,:ni][p][:,p]
        point['physical_banded_eigenvalues']=banded_spectrum(physical_scaled)
        result['points'].append(point)
    result['indicative_probe_wall_s']=time.perf_counter()-started
    return result

def scaling_cases(n):
    d,m=read_input(INPUTS/f'n{n}-f40.npz')
    with np.load(ARCHIVE/f'n{n}-f40-krylov_controle-passage0/bases.npz',allow_pickle=False) as z:
        ii=z['indices_interieur']; b=z['B']
    k=(d[:,ii].T@d[:,ii]).tocsr(); mi=m[ii][:,ii]; ni=len(ii)
    sp=1/np.sqrt(k.diagonal())
    bn=float(np.sqrt(np.sum(b[:,0]**2/mi.diagonal())))
    result=[]
    for hz in [60,80]:
        gamma=(2*np.pi*hz)**2
        c=bmat([[k-gamma*mi,csr_matrix(b)],[csr_matrix(b.T),None]],format='csr')
        choices=[('none',np.ones(ni),1.),('physical_only',sp,1.),
                 ('physical_pow2_only',2.**np.rint(np.log2(sp)),1.),
                 ('inertial',sp,np.sqrt(gamma)/bn),
                 ('inertial_pow2',2.**np.rint(np.log2(sp)),2.**np.rint(np.log2(np.sqrt(gamma)/bn)))]
        for label,pp,bb in choices:
            scale=np.r_[pp,bb]
            a=(diags(scale)@c@diags(scale)).tocsr()
            started=time.perf_counter()
            probe=scalar_ldl(a)
            probe.update(n=n,hz=hz,scaling=label,border_scale=bb,
                         physical_scale_range=[float(pp.min()),float(pp.max())],
                         indicative_ldl_and_reconstruction_s=time.perf_counter()-started)
            result.append(probe)
    return result

def main():
    report=dict(status='SONDE NON CERTIFIEE, binary64; not a comparative campaign',
                date_utc=datetime.now(timezone.utc).isoformat(),
                source_sha256=sha(Path(__file__)),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
                platform=platform.platform(),threads={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS')},cases=[],scalings=[])
    for n in [32,128,512]:
        result=case(n)
        report['cases'].append(result)
        report['scalings'].extend(scaling_cases(n))
        print(json.dumps(dict(n=n,structure=result['structure'],points=[dict(hz=p['hz'],orders=p['orders']) for p in result['points']])),flush=True)
        (OUTPUT/'rapport.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    main()

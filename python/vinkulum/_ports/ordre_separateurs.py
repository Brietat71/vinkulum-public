"""Dissection BFS déterministe du motif physique ; multiplicateurs conservés à la fin.

Heuristique générique, sans géométrie ni coefficient matériel spécifique.
Une permutation ne prouve ni le remplissage ni la stabilité des pivots.
"""
import numpy as np
from .matrices_certificat import matrice_entree
from vinkulum._ports.inverse_selectionnee import _entier_positif


def ordre_separateurs(d,m,*,taille_feuille=8):
    taille_feuille = _entier_positif(taille_feuille,'taille_feuille')
    d,m = matrice_entree(d,'D_i'),matrice_entree(m,'M_ii')
    if m.shape != (d.shape[1],d.shape[1]):
        raise ValueError('dimensions physiques compatibles requises')
    # Produit booléen : couvre le Gram exact même si des termes se compensent.
    pattern = d.astype(bool)
    graph = ((pattern.T@pattern)+m.astype(bool)+m.T.astype(bool)).tocsr()
    graph.setdiag(False);graph.eliminate_zeros()
    adjacency = [set(map(int,graph.indices[graph.indptr[i]:graph.indptr[i+1]])) for i in range(graph.shape[0])]
    def bfs(nodes,start):
        visited={start};layers=[[start]]
        while True:
            frontier=set()
            for v in layers[-1]:
                frontier.update(adjacency[v]&nodes)
            frontier-=visited
            if not frontier:
                return layers,visited
            layers.append(sorted(frontier));visited|=frontier
    order=[];stack=[set(range(graph.shape[0]))]
    while stack:
        nodes=stack.pop()
        if isinstance(nodes,list):
            order.extend(nodes)
        elif len(nodes)<=taille_feuille:
            order.extend(sorted(nodes))
        else:
            layers,visited=bfs(nodes,min(nodes))
            if visited!=nodes:
                stack.extend((nodes-visited,visited))
                continue
            layers,_=bfs(nodes,min(layers[-1]))
            counts=np.cumsum([len(layer) for layer in layers])
            mid=int(np.searchsorted(counts,len(nodes)/2))
            separator=set(layers[mid])
            left=set(v for layer in layers[:mid] for v in layer)
            right=nodes-separator-left
            stack.extend((sorted(separator),right,left))
    return np.array(order,dtype=np.int64)

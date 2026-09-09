# Contraintes redondantes : quotient exact et unicité mécanique

Depuis **0.14.0**, `certifier_initialisation(noyau, redondances=True)`
peut certifier une initialisation avec des contraintes redondantes.
La sortie `vinkulum.quotient.1` est relue par `verifier_certificat` ou
par son script autonome, avec la bibliothèque standard Python seule.

Le contrat porte sur les nombres exportés, après calcul de la géométrie.
Il établit le rang exact des contraintes, leur compatibilité au niveau
des accélérations, l'accélération unique et la réaction généralisée unique.
Il ne déclare pas uniques des multiplicateurs qui ne le sont pas.

## 1. Système complet, base et choix de représentant

Considérons les équations stockées, interprétées exactement :

    M a + Gᵀ λ = f
    G a = c

a contient les accélérations, λ les multiplicateurs. Le noyau a retenu
un ensemble I de r lignes actives. Notons C = G[I,:] et c_I = c[I].
La certification exige un témoin rationnel T tel que :

    G = T C          c = T c_I          T[I,:] = I_r

Ces égalités sont vérifiées **exactement**, sans tolérance de rang ou de
compatibilité. Des coefficients comme 1/3 restent des fractions ; ils ne
sont pas arrondis en binary64.

Le système effectivement résolu par le noyau impose λ_j=0 pour les
lignes inactives. Après permutation, il est la somme directe de :

    K = [ M   Cᵀ ]          [a]   [ f ]
        [ C    0 ]      K · [μ] = [c_I]

et d'un bloc identité sur les multiplicateurs inactifs. Le document
embarque le [certificat linéaire de ce système avec jauge](CERTIFICATION_LINEAIRE.md).
Le vérificateur contrôle sa structure, son second membre et ses
contributions ; une matrice inversible quelconque ne peut pas lui être
substituée sans être détectée.

L'inversibilité est établie par le critère exact `||I−RA||∞<1` de la
0.13.0, pour la matrice avec jauge A réellement exportée. Elle implique
celle de K. La positivité de M n'est pas nécessaire au théorème ci-dessous
et n'est pas revendiquée comme une nouvelle propriété certifiée ici.

## 2. Démonstration

**Rang.** Si les lignes de C étaient dépendantes, il existerait v≠0 avec
Cᵀv=0. Alors K(0,v)=0, contrairement à l'inversibilité de K. C a donc
rang r. Comme G=T C et C est une sous-matrice de G, rang(G)=r exactement.

**Existence.** La solution unique (a*,μ*) de K définit un multiplicateur
λ* par λ*[I]=μ* et λ*[hors I]=0. On a G a*=T C a*=T c_I=c et
Gᵀλ*=Cᵀμ*, donc (a*,λ*) satisfait le système complet.

**Unicité de l'accélération.** Pour toute solution (a,λ) du système
complet, poser μ=Tᵀλ donne M a+Cᵀμ=f et C a=c_I. L'inversibilité de K
impose a=a* et μ=μ*. Ainsi l'accélération est unique, même si λ ne l'est pas.

**Réaction généralisée.** Pour toute solution, Gᵀλ=f−M a*. Cette réaction
est donc unique. Elle désigne le terme d'effort placé à gauche dans
l'équilibre écrit ci-dessus ; avec une convention d'effort appliqué
opposée, le signe change mais les bornes absolues sont les mêmes.

**Multiplicateurs.** Toutes les solutions sont

    λ = λ* + ν,     ν ∈ ker(Gᵀ) = ker(Tᵀ).

En effet, Cᵀ est injective. La dimension de cet espace libre vaut m−r.
La répartition individuelle des réactions de liaison est unique seulement
si r=m. Le représentant choisi par la jauge native n'est pas nécessairement
celui de norme euclidienne minimale ; aucune telle propriété n'est annoncée.

**Erreur.** Le sous-certificat linéaire donne, pour le vecteur calculé
effectivement exporté (â,λ̂), des bornes composantes h :

    |â_i − a*_i| ≤ h_i
    |λ̂_j − λ*_j| ≤ h_(n+j)

Les bornes sur les multiplicateurs concernent le représentant λ* défini
ci-dessus. Elles n'imposent pas que tous les autres représentants soient
proches. Pour la réaction généralisée évaluée exactement sur λ̂ :

    |(Gᵀλ̂ − Gᵀλ*)_i| ≤ Σ_j |G_ji| h_(n+j).

Le générateur arrondit ce dernier majorant vers le haut. Le vérificateur
recalcule la somme en rationnels et refuse une borne insuffisante. Une
évaluation flottante ultérieure de Gᵀλ̂ ajoute ses propres arrondis, qui
ne sont pas inclus dans cette expression.

## 3. Export natif et vérification indépendante

Le chemin natif capture toutes les contributions de G avant filtrage,
le second membre c complet, le masque actif, ainsi que le système et
le vecteur effectivement résolus. Les coefficients actifs utilisés par
le solveur et ceux conservés pour la preuve proviennent de la même
évaluation de G. Les contributions sont additionnées exactement en Python ;
si leur somme n'est pas représentable exactement en binary64, le certificat
est refusé. Le document conserve ces contributions pour relecture.

Le générateur obtient T par élimination de Gauss-Jordan rationnelle sur
les lignes sélectionnées. Le vérificateur ne reprend pas cette élimination :
il contrôle directement les produits T C et T c_I. Il vérifie aussi le
sous-certificat linéaire par le calcul rationnel indépendant de la 0.13.0,
la structure du système avec jauge et les bornes de réaction généralisée.

Le QR flottant du noyau propose le masque. Son seuil de 10⁻¹⁰ ne constitue
pas une preuve de dépendance exacte. Une ligne écartée mais indépendante
en arithmétique exacte fait refuser le certificat, même si le noyau accepte
son initialisation numérique ordinaire. Aucun seuil n'est élargi pour
faire passer cette contre-épreuve.

Le modèle est copié avant évaluation. Les tests confrontent les accélérations
et réactions exportées à celles enregistrées au départ du premier pas de
simulation et vérifient la conservation de l'état lors d'un refus.

## 4. Exemple : cinq contraintes dupliquées

```python
import json
import numpy as np
from vinkulum import Noyau
from vinkulum.certification import certifier_initialisation, verifier_certificat

n = Noyau([0., -9.81, 0.])
n.corps("pendule", 1., (.2*np.eye(3)).ravel().tolist(),
        [0., -1., 0.], v=[.3, 0., 0.], w=[0., 0., .3])
for nom in ("pivot", "copie"):
    n.liaison(nom, None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
n.assemble()
preuve = certifier_initialisation(n, redondances=True, erreur_max=1e-11)
controle = verifier_certificat(preuve)
print(controle["rang_contraintes"], controle["nombre_contraintes"])
print(controle["multiplicateurs_uniques"])
with open("quotient.json", "w") as f:
    json.dump(preuve, f, indent=2)
```

Résultat : **rang 5 sur 10**, multiplicateurs non uniques, accélération
et réaction généralisée uniques. Une [preuve complète est archivée](bancs/certificat-quotient-pendule-0.14.0.json).

```bash
python3 -I -S python/vinkulum/_verification_lineaire.py quotient.json
```

`erreur_max` contrôle la borne numérique globale sur les accélérations
et le représentant des multiplicateurs, dans leurs unités respectives.
Les majorants des réactions généralisées sont fournis séparément dans
`bornes_reaction_generalisee` ; ce paramètre ne leur impose pas un seuil.
La norme globale n'est pas une norme physique invariante au choix d'unités.

## 5. Contre-épreuves et limites

La campagne conserve huit certificats : six pendules redondants sur trois
échelles, un corps libre et une liaison non holonome dupliquée. Les critères
de précision et de rang sont fixés avant exécution dans
`ci/audit_certification_quotient.py`.

Les tests couvrent aussi des dépendances rationnelles non dyadiques,
permutations, grandes variations d'échelle, rang zéro et comparaison à une
solution rationnelle indépendante. Ils altèrent T, G, c, les contributions,
le masque, la jauge et les majorants, puis exigent le refus.

Deux pièges sont explicitement contrôlés :

- Deux lignes `[1,0]` et `[1,2⁻⁸⁰]` sont indépendantes malgré leur proximité.
  Un exemple natif emploie deux points de liaison séparés par ce petit
  bras de levier ; le QR en écarte une ligne, mais le certificat refuse.
- Deux lignes identiques avec seconds membres 0 et 2⁻¹⁰⁷⁴ sont exactement
  incompatibles. Leur écart ne devient pas acceptable par sous-débordement.

L'option reste facultative. La limite de 128 inconnues porte sur le
système exporté **avant retrait des lignes redondantes**, avec un budget
par défaut de 64. Les contacts non lisses, les partitions gelées et les
schémas de redémarrage imposés restent exclus. L'arithmétique rationnelle
peut devenir coûteuse ; aucune qualification des grands systèmes n'est
annoncée. Sans l'option, le contrat de la 0.13.0 est conservé.

Le rang exact des nombres stockés n'est pas nécessairement celui des
fonctions mécaniques idéales. Une dépendance exacte peut être détruite par
l'arrondi des coefficients. Les redondances avec perturbations, l'incertitude
des coefficients, les trajectoires et les événements restent à certifier.

La base de confiance reste celle du certificat linéaire, augmentée de la
fidélité du témoin des contraintes complètes. Les inégalités et identités
sont mathématiquement démontrées et exécutables ; le code n'est pas
formellement vérifié dans un assistant de preuve.

## 6. Position scientifique

[Rump (2013), *Improved componentwise verified error bounds for least squares
problems and underdetermined linear systems*](https://www.tuhh.de/ti3/paper/rump/Ru13a.pdf)
souligne la discontinuité de la pseudo-inverse à rang déficient et établit
des bornes sous une hypothèse de plein rang vérifiée a posteriori. Cette
limite explique pourquoi un seuil numérique ne peut pas certifier ici
une égalité de rang. La présente extension utilise une factorisation exacte
des dépendances, puis un système quotient inversible ; elle ne prétend
implémenter l'algorithme de cet article ni résoudre le problème perturbé.

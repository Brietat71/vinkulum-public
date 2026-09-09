import Mathlib.Data.Rat.Cast.Order

namespace Vinkulum

/-- `puissance` est l’exposant dyadique décalé de 1074. -/
structure Dyadique where
  negatif : Bool
  mantisse : ℕ
  puissance : ℕ
  deriving Repr

def signe (negatif : Bool) (n : ℕ) : ℤ := if negatif then -(n : ℤ) else n

def valeur (d : Dyadique) : ℚ :=
  (signe d.negatif d.mantisse : ℚ) * 2 ^ d.puissance / 2 ^ 1074

def uniteInverse : ℚ := 2 ^ 2148

def coefficient (d : Dyadique) : ℤ :=
  signe d.negatif d.mantisse * 2 ^ (d.puissance + 1074)

def produit (a b : Dyadique) : ℤ :=
  signe (a.negatif != b.negatif) (a.mantisse * b.mantisse) *
    2 ^ (a.puissance + b.puissance)

def residu (termes : List (Dyadique × Dyadique)) (dt : Dyadique) : ℤ :=
  coefficient dt + (termes.map fun p => produit p.1 p.2).sum

def echelle (termes : List (Dyadique × Dyadique)) (dt : Dyadique) : ℕ :=
  (coefficient dt).natAbs + (termes.map fun p => (produit p.1 p.2).natAbs).sum

def residuMathematique (termes : List (Dyadique × Dyadique)) (dt : Dyadique) : ℚ :=
  valeur dt + (termes.map fun p => valeur p.1 * valeur p.2).sum

def echelleMathematique (termes : List (Dyadique × Dyadique)) (dt : Dyadique) : ℚ :=
  |valeur dt| + (termes.map fun p => |valeur p.1 * valeur p.2|).sum

def accepte (termes : List (Dyadique × Dyadique)) (dt tol : Dyadique) : Bool :=
  decide ((residu termes dt).natAbs * 2 ^ 46 ≤
    (coefficient tol).natAbs * 2 ^ 46 + echelle termes dt)

def accumule (termes : List (Dyadique × Dyadique)) (etat : ℤ × ℕ) : ℤ × ℕ :=
  termes.foldl (fun s p => (s.1 + produit p.1 p.2, s.2 + (produit p.1 p.2).natAbs)) etat

abbrev Mot64 := Fin (2 ^ 64)

def champ (w : Mot64) : ℕ := w.val / 2 ^ 52 % 2 ^ 11

def fraction (w : Mot64) : ℕ := w.val % 2 ^ 52

def negatif (w : Mot64) : Bool := decide (2 ^ 63 ≤ w.val)

def decode (w : Mot64) : Option Dyadique :=
  if champ w = 2047 then none else
    some ⟨negatif w, fraction w + (if champ w = 0 then 0 else 2 ^ 52),
      if champ w = 0 then 0 else champ w - 1⟩

/-- Significande IEEE et exposant biaisé, exprimés sans arrondi. -/
def valeurIEEE (w : Mot64) : ℚ :=
  (signe (negatif w) 1 : ℚ) *
    if champ w = 0 then ((fraction w : ℚ) / 2 ^ 52) / 2 ^ 1022
    else (1 + (fraction w : ℚ) / 2 ^ 52) * 2 ^ champ w / 2 ^ 1023

end Vinkulum

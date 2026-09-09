import Mathlib
import VinkulumPreuves.Modele

/-!
Modèle mathématique du garde de vitesse. Les entiers sont non bornés.
La correspondance avec num-bigint/Rust reste une obligation distincte.
-/
set_option maxRecDepth 4096
set_option maxHeartbeats 2000000

namespace Vinkulum

theorem signe_produit (a b : Bool) (m n : ℕ) :
    signe (a != b) (m * n) = signe a m * signe b n := by
  cases a <;> cases b <;> simp [signe]

theorem unite_positive : 0 < uniteInverse := by norm_num [uniteInverse]

theorem coefficient_exact (d : Dyadique) :
    (coefficient d : ℚ) = valeur d * uniteInverse := by
  simp only [coefficient, valeur, uniteInverse, Int.cast_mul, Int.cast_pow, Int.cast_ofNat]
  rw [pow_add]
  norm_num
  ring

theorem produit_exact (a b : Dyadique) :
    (produit a b : ℚ) = (valeur a * valeur b) * uniteInverse := by
  rw [produit, signe_produit]
  simp only [valeur, uniteInverse, Int.cast_mul, Int.cast_pow, Int.cast_ofNat]
  rw [pow_add]
  norm_num
  ring

theorem somme_exacte (termes : List (Dyadique × Dyadique)) :
    (((termes.map fun p => produit p.1 p.2).sum : ℤ) : ℚ) =
      (termes.map fun p => valeur p.1 * valeur p.2).sum * uniteInverse := by
  induction termes with
  | nil => simp
  | cons p ps ih => simp [ih, produit_exact, add_mul]

theorem residu_exact (termes : List (Dyadique × Dyadique)) (dt : Dyadique) :
    (residu termes dt : ℚ) = residuMathematique termes dt * uniteInverse := by
  simp only [residu, residuMathematique, Int.cast_add]
  rw [coefficient_exact, somme_exacte, add_mul]

theorem magnitude_exacte (i : ℤ) (v : ℚ) (h : (i : ℚ) = v * uniteInverse) :
    (i.natAbs : ℚ) = |v| * uniteInverse := by
  rw [Int.cast_natAbs, Int.cast_abs, h, abs_mul, abs_of_pos unite_positive]

theorem echelle_exacte (termes : List (Dyadique × Dyadique)) (dt : Dyadique) :
    (echelle termes dt : ℚ) = echelleMathematique termes dt * uniteInverse := by
  have hc := magnitude_exacte _ _ (coefficient_exact dt)
  have hs : (((termes.map fun p => (produit p.1 p.2).natAbs).sum : ℕ) : ℚ) =
      (termes.map fun p => |valeur p.1 * valeur p.2|).sum * uniteInverse := by
    induction termes with
    | nil => simp
    | cons p ps ih =>
      simp [ih, magnitude_exacte _ _ (produit_exact p.1 p.2), add_mul]
  simpa [echelle, echelleMathematique, add_mul] using congrArg₂ (· + ·) hc hs

theorem garde_exact (termes : List (Dyadique × Dyadique)) (dt tol : Dyadique)
    (htol : 0 ≤ valeur tol) :
    accepte termes dt tol = true ↔
      |residuMathematique termes dt| ≤ valeur tol + echelleMathematique termes dt / 2 ^ 46 := by
  rw [accepte, decide_eq_true_eq]
  have hr := magnitude_exacte _ _ (residu_exact termes dt)
  have ht := magnitude_exacte _ _ (coefficient_exact tol)
  rw [abs_of_nonneg htol] at ht
  have he := echelle_exacte termes dt
  have cast : ((residu termes dt).natAbs * 2 ^ 46 ≤
      (coefficient tol).natAbs * 2 ^ 46 + echelle termes dt) ↔
      ((residu termes dt).natAbs : ℚ) * 2 ^ 46 ≤
        ((coefficient tol).natAbs : ℚ) * 2 ^ 46 + (echelle termes dt : ℚ) := by
    exact_mod_cast Iff.rfl
  rw [cast, hr, ht, he]
  have qpos : (0 : ℚ) < 2 ^ 46 := by positivity
  calc
    |residuMathematique termes dt| * uniteInverse * 2 ^ 46 ≤
        valeur tol * uniteInverse * 2 ^ 46 + echelleMathematique termes dt * uniteInverse
      ↔ (|residuMathematique termes dt| * 2 ^ 46) * uniteInverse ≤
          (valeur tol * 2 ^ 46 + echelleMathematique termes dt) * uniteInverse := by ring_nf
    _ ↔ |residuMathematique termes dt| * 2 ^ 46 ≤
          valeur tol * 2 ^ 46 + echelleMathematique termes dt :=
      mul_le_mul_iff_of_pos_right unite_positive
    _ ↔ |residuMathematique termes dt| ≤
          valeur tol + echelleMathematique termes dt / 2 ^ 46 := by
      rw [← le_div_iff₀ qpos]
      field_simp

/-- Accumulation séquentielle correspondant à la boucle du garde. -/
theorem accumule_exact (termes : List (Dyadique × Dyadique)) (r : ℤ) (e : ℕ) :
    accumule termes (r, e) =
      (r + (termes.map fun p => produit p.1 p.2).sum,
       e + (termes.map fun p => (produit p.1 p.2).natAbs).sum) := by
  induction termes generalizing r e with
  | nil => simp [accumule]
  | cons p ps ih =>
    simp only [accumule, List.foldl_cons]
    change accumule ps (r + produit p.1 p.2, e + (produit p.1 p.2).natAbs) = _
    rw [ih]
    simp [add_assoc]

theorem boucle_exacte (termes : List (Dyadique × Dyadique)) (dt : Dyadique) :
    accumule termes (coefficient dt, (coefficient dt).natAbs) =
      (residu termes dt, echelle termes dt) := by
  exact accumule_exact termes _ _

theorem decalage_exact (r t e : ℕ) :
    (r <<< 46 ≤ (t <<< 46) + e) ↔ r * 2 ^ 46 ≤ t * 2 ^ 46 + e := by
  simp only [Nat.shiftLeft_eq]

end Vinkulum

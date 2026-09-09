import VinkulumPreuves.GardeVitesse

set_option maxRecDepth 4096
set_option maxHeartbeats 2000000
namespace Vinkulum

theorem fraction_bornee (w : Mot64) : fraction w < 2 ^ 52 :=
  Nat.mod_lt _ (by positivity)
theorem champ_borne (w : Mot64) : champ w < 2048 :=
  Nat.mod_lt _ (by norm_num)
theorem bit_implicite (w : Mot64) :
    fraction w ||| 2 ^ 52 = fraction w + 2 ^ 52 := by
  have h := Nat.two_pow_add_eq_or_of_lt (fraction_bornee w) 1
  simpa [Nat.or_comm, Nat.add_comm] using h.symm

theorem signe_unite (b : Bool) (n : ℕ) :
    (signe b n : ℚ) = (signe b 1 : ℚ) * n := by
  cases b <;> simp [signe]

theorem decode_exact (w : Mot64) (d : Dyadique) (h : decode w = some d) :
    valeur d = valeurIEEE w := by
  unfold decode at h
  split at h
  · contradiction
  · cases h
    simp only [valeur, valeurIEEE]
    split_ifs with hz
    · rw [signe_unite]
      simp only [Nat.add_zero, pow_zero, mul_one]
      norm_num
      ring
    · have he : champ w - 1 + 1 = champ w := by omega
      have hp : (2 : ℚ) ^ champ w = 2 ^ (champ w - 1) * 2 := by
        conv_lhs => rw [← he]
        rw [pow_add]
        norm_num
      rw [signe_unite]
      simp only [Nat.cast_add, Nat.cast_pow, Nat.cast_ofNat]
      rw [hp]
      norm_num
      ring

theorem decode_bornes (w : Mot64) (d : Dyadique) (h : decode w = some d) :
    d.mantisse < 2 ^ 53 ∧ d.puissance ≤ 2045 := by
  have hf := fraction_bornee w
  have he := champ_borne w
  unfold decode at h
  split at h
  · contradiction
  · cases h
    simp only
    split_ifs <;> norm_num at * <;> omega

theorem produit_u128 (w v : Mot64) (a b : Dyadique)
    (ha : decode w = some a) (hb : decode v = some b) :
    a.mantisse * b.mantisse < 2 ^ 106 ∧ a.mantisse * b.mantisse < 2 ^ 128 := by
  have h := Nat.mul_lt_mul_of_lt_of_lt (decode_bornes w a ha).1 (decode_bornes v b hb).1
  norm_num at h ⊢
  constructor <;> omega

theorem exposants_rust (w v : Mot64) (a b : Dyadique)
    (ha : decode w = some a) (hb : decode v = some b) :
    -2148 ≤ ((a.puissance : ℤ) - 1074) + ((b.puissance : ℤ) - 1074) ∧
    ((a.puissance : ℤ) - 1074) + ((b.puissance : ℤ) - 1074) ≤ 1942 ∧
    ((a.puissance : ℤ) - 1074) + ((b.puissance : ℤ) - 1074) + 2148 =
      (a.puissance + b.puissance : ℕ) := by
  have ha' := (decode_bornes w a ha).2
  have hb' := (decode_bornes v b hb).2
  omega

theorem tolerance_admise (d : Dyadique)
    (h : (d.negatif && d.mantisse != 0) = false) : 0 ≤ valeur d := by
  cases d with
  | mk neg m p =>
    cases neg
    · simp only [valeur, signe, Bool.false_eq_true, if_false, Int.cast_natCast]
      positivity
    · simp_all [valeur, signe]

theorem fraction_masque (w : Mot64) :
    w.val &&& (2 ^ 52 - 1) = fraction w := by
  exact Nat.and_two_pow_sub_one_eq_mod _ _

theorem champ_masque (w : Mot64) :
    (w.val >>> 52) &&& (2 ^ 11 - 1) = champ w := by
  simp only [Nat.and_two_pow_sub_one_eq_mod, Nat.shiftRight_eq_div_pow, champ]

theorem signe_decalage (w : Mot64) :
    (w.val >>> 63 != 0) = negatif w := by
  have h : w.val / 2 ^ 63 ≠ 0 ↔ 2 ^ 63 ≤ w.val := by
    have hz := Nat.div_eq_zero_iff (a := w.val) (b := 2 ^ 63)
    omega
  apply Bool.eq_iff_iff.mpr
  simpa only [Nat.shiftRight_eq_div_pow, negatif, bne_iff_ne, decide_eq_true_eq] using h

theorem refus_non_fini (w : Mot64) : decode w = none ↔ champ w = 2047 := by
  simp only [decode]
  split <;> simp_all

end Vinkulum

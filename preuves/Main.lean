import VinkulumPreuves.Modele

open Vinkulum

def lireMot (s : String) : Except String Dyadique := do
  let n ← match s.toNat? with | some n => pure n | none => throw "entier"
  if h : n < 2 ^ 64 then
    match decode ⟨n, h⟩ with | some d => pure d | none => throw "non_fini"
  else throw "mot_hors_plage"

def paires : List Dyadique → Except String (List (Dyadique × Dyadique))
  | [] => .ok []
  | a :: b :: xs => (paires xs).map ((a, b) :: ·)
  | _ => .error "dimensions"

def calcule (line : String) : Except String String := do
  let mots := line.trim.splitOn " "
  let ds ← mots.mapM lireMot
  match ds with
  | tol :: dt :: xs =>
    if tol.negatif && tol.mantisse != 0 then throw "tolerance_negative"
    let ts ← paires xs
    let (r, e) := accumule ts (coefficient dt, (coefficient dt).natAbs)
    let t := (coefficient tol).natAbs
    let ok := decide (r.natAbs <<< 46 ≤ (t <<< 46) + e)
    return s!"{ok} {r} {e} {t}"
  | _ => throw "dimensions"

def main : IO Unit := do
  let stdin ← IO.getStdin
  let stdout ← IO.getStdout
  repeat
    let line ← stdin.getLine
    if line.isEmpty then break
    match calcule line with
    | .ok s => stdout.putStrLn s
    | .error e => stdout.putStrLn s!"ERREUR {e}"

//! Prototype statique indépendant : poutre mixte hybride d'ordres 1 et 2.
//! Dérivation propre d'après Humer–Steinbrecher–Pechstein (2026), équation 49.
//! https://arxiv.org/html/2605.04573v1 — aucun code de solveur tiers utilisé.
//! Ce programme n'ajoute pas de formulation à l'API distribuée.

#[allow(dead_code)]
#[path = "../src/ad.rs"]
mod ad;
#[path = "poutre_mixte/element.rs"]
mod element;
#[path = "poutre_mixte/systeme.rs"]
mod systeme;
#[cfg(test)]
#[path = "poutre_mixte/tests.rs"]
mod tests;

use nalgebra::{Matrix3, Vector3};
type V3 = Vector3<f64>;
type M3 = Matrix3<f64>;

fn main() {
    // Les jets locaux emboîtés sont volumineux ; pile bornée explicitement
    // pour le prototype, sans changer la pile du noyau Python.
    let result = std::thread::Builder::new()
        .stack_size(16 * 1024 * 1024)
        .spawn(systeme::main)
        .expect("thread de calcul")
        .join()
        .expect("calcul sans panique");
    if let Err(error) = result {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

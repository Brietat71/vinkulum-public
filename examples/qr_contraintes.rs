//! Diagnostic du QR creux installé sur des contraintes presque dépendantes.
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::linalg::qr::{factorize_symbolic_qr, QrSymbolicParams};
use faer::sparse::linalg::SupernodalThreshold;
use faer::sparse::{SparseColMat, Triplet};

fn probe(delta: f64, supernodal: bool) {
    let entries = [
        Triplet {
            row: 0,
            col: 0,
            val: 1.0,
        },
        Triplet {
            row: 0,
            col: 1,
            val: 1.0,
        },
        Triplet {
            row: 1,
            col: 1,
            val: delta,
        },
    ];
    let a = SparseColMat::<usize, f64>::try_new_from_triplets(8, 2, &entries).unwrap();
    let sym = factorize_symbolic_qr(
        a.symbolic(),
        QrSymbolicParams {
            supernodal_flop_ratio_threshold: if supernodal {
                SupernodalThreshold::FORCE_SUPERNODAL
            } else {
                SupernodalThreshold::FORCE_SIMPLICIAL
            },
            ..Default::default()
        },
    )
    .unwrap();
    let mut indices = vec![0; sym.len_idx()];
    let mut values = vec![0.0; sym.len_val()];
    let mut mem =
        MemBuffer::new(sym.factorize_numeric_qr_scratch::<f64>(faer::Par::Seq, Default::default()));
    let qr = sym.factorize_numeric_qr(
        &mut indices,
        &mut values,
        a.as_ref(),
        faer::Par::Seq,
        MemStack::new(&mut mem),
        Default::default(),
    );
    let mut rhs = faer::Mat::zeros(8, 1);
    rhs[(0, 0)] = 2.0;
    rhs[(1, 0)] = 1.0;
    rhs[(2, 0)] = 3.0;
    let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<f64>(1, faer::Par::Seq));
    qr.solve_in_place_with_conj(
        faer::Conj::No,
        rhs.as_mut(),
        faer::Par::Seq,
        MemStack::new(&mut mem),
    );
    println!(
        "delta={delta:e} supernodal={supernodal} lambda={:?} r={:?}",
        [rhs[(0, 0)], rhs[(1, 0)]],
        [
            2.0 - rhs[(0, 0)] - rhs[(1, 0)],
            1.0 - delta * rhs[(1, 0)],
            3.0
        ]
    );
}

fn main() {
    faer::set_global_parallelism(faer::Par::Seq);
    for delta in [0.0, 1e-12, 1e-9, 1e-6, 1.0] {
        for supernodal in [false, true] {
            let result = std::panic::catch_unwind(|| probe(delta, supernodal));
            if result.is_err() {
                println!("PANIC delta={delta:e} supernodal={supernodal}");
            }
        }
    }
}

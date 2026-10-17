use std::sync::atomic::AtomicBool;

use graynfs_core::algebraic::{QuadraticElement, QuadraticRing};
use graynfs_core::base::FactorBase;
use graynfs_core::integer::Integer;
use graynfs_core::linear_algebra::{find_dependencies, matrix_rank};
use graynfs_core::options::FactorOptions;
use graynfs_core::polynomial::select_polynomial;
use graynfs_core::relations::collect_relations;
use graynfs_core::square_root::extract_factor;

#[test]
fn inspect_quadratic_relations() {
    let n = graynfs_core::parse_decimal("1000036000099").expect("valid integer");
    let options = FactorOptions {
        degree: Some(2),
        sieve_a_bound: 700,
        sieve_b_bound: 120,
        rational_base_bound: 2_000,
        algebraic_base_bound: 2_000,
        max_duration: None,
        ..FactorOptions::default()
    };
    let polynomial = select_polynomial(&n, &options).expect("polynomial");
    let base = FactorBase::build(&polynomial, &n, 2_000, 2_000).expect("base");
    let mut progress = |_: &graynfs_core::ProgressEvent| Ok(());
    let relations = collect_relations(
        &n,
        &polynomial,
        &base,
        &options,
        &AtomicBool::new(false),
        None,
        &mut progress,
    )
    .expect("relations");
    let rank = matrix_rank(&relations, base.column_count());
    let dependencies = find_dependencies(&relations, base.column_count(), 128);
    println!(
        "base={} algebraic={} relations={} rank={} dependencies={}",
        base.rational_primes.len(),
        base.algebraic_primes.len(),
        relations.len(),
        rank,
        dependencies.len()
    );
    for (dependency_index, dependency) in dependencies.iter().enumerate() {
        let ring = QuadraticRing::new(polynomial.constant());
        let mut rational_product = graynfs_core::Integer::from(1);
        let mut algebraic_product = QuadraticElement {
            x: graynfs_core::Integer::from(1),
            y: graynfs_core::Integer::from(0),
        };
        let mut combined_mask = vec![0_u64; base.column_count().div_ceil(64)];
        for index in dependency {
            let relation = &relations[*index];
            let set_bits: Vec<usize> = relation
                .mask
                .iter()
                .enumerate()
                .flat_map(|(word_index, word)| {
                    (0..64).filter_map(move |bit| {
                        (word & (1_u64 << bit) != 0).then_some(word_index * 64 + bit)
                    })
                })
                .collect();
            println!(
                "relation index={} a={} b={} rational={} algebraic=({}, {}) mask={} bits={:?}",
                index,
                relation.a,
                relation.b,
                relation.rational_value,
                relation.algebraic_value.x,
                relation.algebraic_value.y,
                relation
                    .mask
                    .iter()
                    .map(|word| word.count_ones())
                    .sum::<u32>(),
                set_bits
            );
            rational_product *= &relation.rational_value;
            algebraic_product = ring.multiply(&algebraic_product, &relation.algebraic_value);
            for (word_index, word) in relation.mask.iter().enumerate() {
                combined_mask[word_index] ^= word;
            }
        }
        let algebraic_root = ring.square_root(&algebraic_product);
        let rational_root = graynfs_core::integer_sqrt(&rational_product);
        let norm = ring.norm(&algebraic_product);
        let norm_root = graynfs_core::integer_sqrt(&norm);
        println!(
            "dependency={} size={} indices={:?} mask={} rational={} rational_square={} norm_square={} algebraic_root={}",
            dependency_index,
            dependency.len(),
            dependency,
            combined_mask
                .iter()
                .map(|word| word.count_ones())
                .sum::<u32>(),
            rational_product,
            &rational_root * &rational_root == rational_product,
            &norm_root * &norm_root == norm,
            algebraic_root.is_some()
        );
        if let Some(root) = algebraic_root {
            println!(
                "product x={} y={} norm={}",
                algebraic_product.x, algebraic_product.y, norm
            );
            println!("root u={} v={}", root.x, root.y);
            for sign_u in [Integer::from(1), Integer::from(-1)] {
                for sign_v in [Integer::from(1), Integer::from(-1)] {
                    let candidate = QuadraticElement {
                        x: &sign_u * &root.x,
                        y: &sign_v * &root.y,
                    };
                    let candidate_square = ring.multiply(&candidate, &candidate);
                    let evaluation = ring.evaluate(&candidate, &polynomial.m);
                    let x = graynfs_core::integer_abs(&(&evaluation % &n));
                    let y = rational_root.clone();
                    let congruence = (&x * &x - &y * &y) % &n;
                    let gcd_minus = graynfs_core::integer_gcd(&(&x - &y), &n);
                    let gcd_plus = graynfs_core::integer_gcd(&(&x + &y), &n);
                    println!(
                        "candidate u={} v={} square={} evaluation={} x={} y={} congruence={} gcd-={} gcd+={}",
                        candidate.x,
                        candidate.y,
                        candidate_square.x == algebraic_product.x
                            && candidate_square.y == algebraic_product.y,
                        evaluation,
                        x,
                        y,
                        congruence,
                        gcd_minus,
                        gcd_plus
                    );
                }
            }
            println!(
                "extract={:?}",
                extract_factor(&n, &ring, &polynomial.m, &relations, dependency)
            );
        }
        if dependency_index >= 9 {
            break;
        }
    }
}

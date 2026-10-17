use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use graynfs_core::algebraic::QuadraticRing;
use graynfs_core::base::FactorBase;
use graynfs_core::linear_algebra::{find_dependencies, matrix_rank};
use graynfs_core::options::FactorOptions;
use graynfs_core::polynomial::{PolynomialPair, select_polynomial};
use graynfs_core::relations::collect_relations;
use graynfs_core::square_root::extract_factor;

fn attempt(number: &graynfs_core::Integer, options: &FactorOptions, polynomial: &PolynomialPair) {
    let base = FactorBase::build(
        polynomial,
        number,
        options.rational_base_bound,
        options.algebraic_base_bound,
    )
    .expect("base");
    let start = Instant::now();
    let mut progress = |_: &graynfs_core::ProgressEvent| Ok(());
    let relations = collect_relations(
        number,
        polynomial,
        &base,
        options,
        &AtomicBool::new(false),
        None,
        &mut progress,
    )
    .expect("relations");
    let rank = matrix_rank(&relations, base.column_count(), 1);
    println!(
        "attempt base={} columns={} relations={} rank={} elapsed={}ms",
        base.rational_primes.len(),
        base.column_count(),
        relations.len(),
        rank,
        start.elapsed().as_millis()
    );
    if relations.len() <= rank {
        return;
    }
    let dependencies = find_dependencies(&relations, base.column_count(), 128, 1);
    println!("dependencies={}", dependencies.len());
    let ring = QuadraticRing::new(polynomial.constant());
    for (di, dependency) in dependencies.iter().enumerate() {
        let mut rational = graynfs_core::Integer::from(1);
        let mut algebraic = graynfs_core::algebraic::QuadraticElement {
            x: graynfs_core::Integer::from(1),
            y: graynfs_core::Integer::from(0),
        };
        for index in dependency {
            let relation = &relations[*index];
            rational *= &relation.rational_value;
            algebraic = ring.multiply(&algebraic, &relation.algebraic_value);
        }
        let norm = ring.norm(&algebraic);
        let norm_root = graynfs_core::integer_sqrt(&norm);
        if &norm_root * &norm_root != norm {
            continue;
        }
        let rational_root = graynfs_core::integer_sqrt(&rational);
        let Some(algebraic_root) = ring.square_root(&algebraic) else {
            println!(
                "dep={} size={} norm_is_square=true square_root=None",
                di,
                dependency.len()
            );
            continue;
        };
        let mut evaluation_product = graynfs_core::Integer::from(1);
        for index in dependency {
            let relation = &relations[*index];
            evaluation_product *=
                &relation.algebraic_value.x - &relation.algebraic_value.y * &polynomial.m;
        }
        println!(
            "m={} alpha_x_digits={} alpha_y_digits={} R_digits={} eval_product_digits={} R==eval_product:{}",
            polynomial.m,
            algebraic.x.to_string().len(),
            algebraic.y.to_string().len(),
            rational.to_string().len(),
            evaluation_product.to_string().len(),
            rational == evaluation_product
        );
        println!(
            "SUCCESS dep={} size={} rational_bits={} norm_bits={} rational_is_square={}",
            di,
            dependency.len(),
            rational.to_string().len(),
            norm.to_string().len(),
            &rational_root * &rational_root == rational
        );
        for sign_u in [1_i64, -1] {
            for sign_v in [1_i64, -1] {
                let candidate = graynfs_core::algebraic::QuadraticElement {
                    x: Integer::from(sign_u) * &algebraic_root.x,
                    y: Integer::from(sign_v) * &algebraic_root.y,
                };
                let square = ring.multiply(&candidate, &candidate);
                let is_root = square.x == algebraic.x && square.y == algebraic.y;
                let evaluation = ring.evaluate(&candidate, &polynomial.m);
                println!(
                    "  variant u={} v={} is_root={} eval_digits={} eval2_minus_R_digits={}",
                    candidate.x.to_string().len(),
                    candidate.y.to_string().len(),
                    is_root,
                    evaluation.to_string().len(),
                    graynfs_core::integer_abs(&(&evaluation * &evaluation - &rational))
                        .to_string()
                        .len()
                );
                let x = &evaluation % number;
                let y = &rational_root % number;
                let congruence = &x * &x - &y * &y;
                println!(
                    "  variant u={} v={} is_root={} eval_bits={} x_bits={} y_bits={} diff_is_zero={} diff_mod_n={} gcd-={} gcd+={}",
                    candidate.x.to_string().len(),
                    candidate.y.to_string().len(),
                    is_root,
                    evaluation.to_string().len(),
                    x.to_string().len(),
                    y.to_string().len(),
                    congruence == graynfs_core::Integer::from(0),
                    &congruence % number == graynfs_core::Integer::from(0),
                    graynfs_core::integer_gcd(&(&x - &y), number),
                    graynfs_core::integer_gcd(&(&x + &y), number)
                );
            }
        }
        println!(
            "  extract={:?}",
            extract_factor(number, &ring, &polynomial.m, &relations, dependency)
        );
        return;
    }
}

use graynfs_core::Integer;

#[test]
#[ignore = "manual quadratic diagnostic"]
fn scratch_quadratic_pipeline() {
    let number = graynfs_core::parse_decimal("1000036000099").expect("valid integer");
    let mut options = FactorOptions {
        degree: Some(2),
        rational_base_bound: 5_000,
        algebraic_base_bound: 5_000,
        sieve_a_bound: 2_500,
        sieve_b_bound: 300,
        max_duration: Some(Duration::from_secs(120)),
        allow_fallback: false,
        ..FactorOptions::default()
    };
    let polynomial = select_polynomial(&number, &options).expect("polynomial");
    for _ in 0..3 {
        attempt(&number, &options, &polynomial);
        options.rational_base_bound = (options.rational_base_bound * 2).min(20_000);
        options.algebraic_base_bound = (options.algebraic_base_bound * 2).min(20_000);
        options.sieve_a_bound = (options.sieve_a_bound * 2).min(8_000);
        options.sieve_b_bound = (options.sieve_b_bound * 2).min(1_000);
    }
}

use crate::relations::Relation;

pub trait Masked {
    fn mask(&self) -> &[u64];
}

impl Masked for Relation {
    fn mask(&self) -> &[u64] {
        &self.mask
    }
}

fn join_candidate_worker(
    handle: std::thread::ScopedJoinHandle<'_, Vec<(usize, usize, Vec<u64>)>>,
) -> Vec<(usize, usize, Vec<u64>)> {
    handle.join().unwrap_or_default()
}

pub fn rank_and_dependencies<M: Masked>(
    relations: &[M],
    column_count: usize,
    limit: usize,
    min_weight: usize,
    workers: usize,
) -> (usize, Vec<Vec<usize>>) {
    if relations.is_empty() || column_count == 0 {
        return (0, Vec::new());
    }
    let mask_words = column_count.div_ceil(64);
    let relation_words = relations.len().div_ceil(64);
    let last_bits = column_count % 64;
    let last_mask = if last_bits == 0 {
        u64::MAX
    } else {
        (1_u64 << last_bits) - 1
    };
    let mut basis_slab = vec![0_u64; column_count * mask_words];
    let mut basis_filled = vec![false; column_count];
    let mut pivot_row = vec![0_usize; column_count];
    let mut order: Vec<usize> = Vec::new();
    let mut used: Vec<Vec<usize>> = Vec::with_capacity(relations.len());
    let mut is_candidate = vec![false; relations.len()];
    let mut rank = 0_usize;
    let mut scratch = vec![0_u64; mask_words];
    for (relation_index, relation) in relations.iter().enumerate() {
        let mask = relation.mask();
        scratch[..mask_words].copy_from_slice(&mask[..mask_words.min(mask.len())]);
        scratch[mask_words - 1] &= last_mask;
        let mut word_index = 0_usize;
        let mut stored = false;
        let mut used_here = Vec::new();
        while word_index < mask_words {
            let mut word = scratch[word_index];
            while word != 0 {
                let column = word_index * 64 + word.trailing_zeros() as usize;
                if column >= column_count {
                    break;
                }
                if basis_filled[column] {
                    let pivot = &basis_slab[column * mask_words..(column + 1) * mask_words];
                    used_here.push(column);
                    for (left, right) in scratch[word_index..]
                        .iter_mut()
                        .zip(pivot[word_index..].iter())
                    {
                        *left ^= *right;
                    }
                    word = scratch[word_index];
                } else {
                    rank += 1;
                    basis_slab[column * mask_words..(column + 1) * mask_words]
                        .copy_from_slice(&scratch);
                    basis_filled[column] = true;
                    pivot_row[column] = relation_index;
                    order.push(column);
                    stored = true;
                    word_index = mask_words;
                    break;
                }
            }
            word_index += 1;
        }
        if stored {
            used.push(used_here);
            continue;
        }
        if limit == 0 {
            used.push(Vec::new());
            continue;
        }
        if scratch.iter().any(|word| *word != 0) {
            used.push(Vec::new());
            continue;
        }
        used.push(used_here);
        is_candidate[relation_index] = true;
    }
    if limit == 0 {
        return (rank, Vec::new());
    }
    let mut combo_slab = vec![0_u64; column_count * relation_words];
    let mut combo_filled = vec![false; column_count];
    for &column in &order {
        let row = pivot_row[column];
        let mut combo = vec![0_u64; relation_words];
        combo[row / 64] |= 1_u64 << (row % 64);
        for &used_column in &used[row] {
            if combo_filled[used_column] {
                let previous =
                    &combo_slab[used_column * relation_words..(used_column + 1) * relation_words];
                for (left, right) in combo.iter_mut().zip(previous.iter()) {
                    *left ^= *right;
                }
            }
        }
        combo_slab[column * relation_words..(column + 1) * relation_words].copy_from_slice(&combo);
        combo_filled[column] = true;
    }
    let candidate_rows: Vec<usize> = (0..used.len()).filter(|row| is_candidate[*row]).collect();
    let mut candidates: Vec<(usize, usize, Vec<u64>)> = Vec::new();
    let workers = workers.max(1);
    if workers == 1 || candidate_rows.len() < 64 {
        for &relation_index in &candidate_rows {
            let mut combo = vec![0_u64; relation_words];
            combo[relation_index / 64] |= 1_u64 << (relation_index % 64);
            for &used_column in &used[relation_index] {
                if combo_filled[used_column] {
                    let previous = &combo_slab
                        [used_column * relation_words..(used_column + 1) * relation_words];
                    for (left, right) in combo.iter_mut().zip(previous.iter()) {
                        *left ^= *right;
                    }
                }
            }
            let weight = combo.iter().map(|word| word.count_ones() as usize).sum();
            if weight >= min_weight {
                candidates.push((weight, relation_index, combo));
            }
        }
    } else {
        let chunk = candidate_rows.len().div_ceil(workers);
        let used = &used;
        let combo_slab = &combo_slab;
        let combo_filled = &combo_filled;
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for rows in candidate_rows.chunks(chunk) {
                handles.push(scope.spawn(move || {
                    let mut local: Vec<(usize, usize, Vec<u64>)> = Vec::new();
                    let mut combo = vec![0_u64; relation_words];
                    for &relation_index in rows {
                        combo.fill(0);
                        combo[relation_index / 64] |= 1_u64 << (relation_index % 64);
                        for &used_column in &used[relation_index] {
                            if combo_filled[used_column] {
                                let previous = &combo_slab[used_column * relation_words
                                    ..(used_column + 1) * relation_words];
                                for (left, right) in combo.iter_mut().zip(previous.iter()) {
                                    *left ^= *right;
                                }
                            }
                        }
                        let weight = combo.iter().map(|word| word.count_ones() as usize).sum();
                        if weight >= min_weight {
                            local.push((weight, relation_index, combo.clone()));
                        }
                    }
                    local
                }));
            }
            for handle in handles {
                candidates.extend(join_candidate_worker(handle));
            }
        });
    }
    if candidates.len() > limit {
        candidates.select_nth_unstable_by_key(limit, |(weight, sequence, _)| (*weight, *sequence));
        candidates.truncate(limit);
    }
    candidates.sort_by_key(|(weight, sequence, _)| (*weight, *sequence));
    let mut dependencies = Vec::with_capacity(candidates.len());
    for (_, _, combination) in &candidates {
        let mut dependency = Vec::new();
        for (word_index, word) in combination.iter().enumerate() {
            let mut rest = *word;
            while rest != 0 {
                let bit = rest.trailing_zeros() as usize;
                let index = word_index * 64 + bit;
                if index < relations.len() {
                    dependency.push(index);
                }
                rest &= rest - 1;
            }
        }
        dependencies.push(dependency);
    }
    (rank, dependencies)
}

pub fn find_dependencies(
    relations: &[Relation],
    column_count: usize,
    limit: usize,
    workers: usize,
) -> Vec<Vec<usize>> {
    rank_and_dependencies(relations, column_count, limit, 1, workers).1
}

pub fn matrix_rank(relations: &[Relation], column_count: usize, workers: usize) -> usize {
    rank_and_dependencies(relations, column_count, 0, usize::MAX, workers).0
}

pub(crate) fn order_dependencies_by_cost(dependencies: &[Vec<usize>]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..dependencies.len()).collect();
    order.sort_by_key(|&index| dependencies[index].len());
    order
}

#[cfg(test)]
mod tests {
    use super::{Masked, join_candidate_worker, order_dependencies_by_cost, rank_and_dependencies};

    struct TestRelation {
        mask: Vec<u64>,
    }

    impl Masked for TestRelation {
        fn mask(&self) -> &[u64] {
            &self.mask
        }
    }

    fn duplicate_relations(count: usize, word: u64) -> Vec<TestRelation> {
        (0..count)
            .map(|_| TestRelation { mask: vec![word] })
            .collect()
    }

    #[test]
    fn panicking_candidate_worker_falls_back_to_empty() {
        std::thread::scope(|scope| {
            let handle = scope.spawn(|| -> Vec<(usize, usize, Vec<u64>)> {
                panic!("injected candidate worker panic");
            });
            let started = std::time::Instant::now();
            let recovered = join_candidate_worker(handle);
            assert!(recovered.is_empty());
            assert_eq!(started.elapsed().as_secs(), 0);
        });
    }

    #[test]
    fn cost_order_sorts_by_length_stably() {
        let dependencies = vec![vec![0, 1, 2], vec![3], vec![4, 5], vec![]];
        assert_eq!(order_dependencies_by_cost(&dependencies), vec![3, 1, 2, 0]);
        let empty: Vec<Vec<usize>> = Vec::new();
        assert!(order_dependencies_by_cost(&empty).is_empty());
    }

    #[test]
    fn single_worker_matches_parallel_workers_on_large_candidate_set() {
        let relations = duplicate_relations(128, 0b0011);
        let (single_rank, single) = rank_and_dependencies(&relations, 4, 128, 1, 1);
        let (parallel_rank, parallel) = rank_and_dependencies(&relations, 4, 128, 1, 8);
        assert_eq!(single_rank, parallel_rank);
        assert_eq!(single, parallel);
        assert!(!single.is_empty());
        for dependency in &single {
            let mut word = 0_u64;
            for &index in dependency {
                word ^= relations[index].mask[0];
            }
            assert_eq!(word, 0);
        }
    }

    #[test]
    fn zero_workers_falls_back_to_single_worker() {
        let relations = duplicate_relations(96, 0b0110);
        let (expected_rank, expected) = rank_and_dependencies(&relations, 4, 128, 1, 1);
        let (actual_rank, actual) = rank_and_dependencies(&relations, 4, 128, 1, 0);
        assert_eq!(expected_rank, actual_rank);
        assert_eq!(expected, actual);
    }
}

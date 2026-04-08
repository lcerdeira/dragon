/// Fenwick tree (Binary Indexed Tree) supporting prefix maximum queries.
/// Used in the colinear chaining algorithm for O(h log h) DP.

pub struct FenwickMax {
    tree: Vec<i64>,
    n: usize,
}

impl FenwickMax {
    /// Create a new Fenwick tree of size n, initialized to i64::MIN.
    pub fn new(n: usize) -> Self {
        Self {
            tree: vec![i64::MIN; n + 1],
            n,
        }
    }

    /// Update position i with value val (keeps the maximum).
    pub fn update(&mut self, mut i: usize, val: i64) {
        i += 1; // 1-indexed
        while i <= self.n {
            self.tree[i] = self.tree[i].max(val);
            i += i & i.wrapping_neg();
        }
    }

    /// Query the maximum value in the prefix [0, i].
    pub fn prefix_max(&self, mut i: usize) -> i64 {
        i += 1; // 1-indexed
        let mut result = i64::MIN;
        while i > 0 {
            result = result.max(self.tree[i]);
            i -= i & i.wrapping_neg();
        }
        result
    }

    /// Reset all values to i64::MIN.
    pub fn reset(&mut self) {
        self.tree.fill(i64::MIN);
    }
}

/// Standard Fenwick tree for prefix sums (used for vote counting).
pub struct FenwickSum {
    tree: Vec<i64>,
    n: usize,
}

impl FenwickSum {
    pub fn new(n: usize) -> Self {
        Self {
            tree: vec![0; n + 1],
            n,
        }
    }

    pub fn update(&mut self, mut i: usize, delta: i64) {
        i += 1;
        while i <= self.n {
            self.tree[i] += delta;
            i += i & i.wrapping_neg();
        }
    }

    pub fn prefix_sum(&self, mut i: usize) -> i64 {
        i += 1;
        let mut result = 0i64;
        while i > 0 {
            result += self.tree[i];
            i -= i & i.wrapping_neg();
        }
        result
    }

    pub fn range_sum(&self, l: usize, r: usize) -> i64 {
        if l == 0 {
            self.prefix_sum(r)
        } else {
            self.prefix_sum(r) - self.prefix_sum(l - 1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fenwick_max() {
        let mut fw = FenwickMax::new(10);
        fw.update(3, 10);
        fw.update(5, 20);
        fw.update(1, 5);

        assert_eq!(fw.prefix_max(0), i64::MIN);
        assert_eq!(fw.prefix_max(1), 5);
        assert_eq!(fw.prefix_max(3), 10);
        assert_eq!(fw.prefix_max(5), 20);
        assert_eq!(fw.prefix_max(9), 20);
    }

    #[test]
    fn test_fenwick_max_overwrite() {
        let mut fw = FenwickMax::new(5);
        fw.update(2, 10);
        fw.update(2, 30);
        assert_eq!(fw.prefix_max(2), 30);
    }

    #[test]
    fn test_fenwick_sum() {
        let mut fw = FenwickSum::new(10);
        fw.update(0, 1);
        fw.update(1, 2);
        fw.update(2, 3);
        assert_eq!(fw.prefix_sum(2), 6);
        assert_eq!(fw.range_sum(1, 2), 5);
    }
}

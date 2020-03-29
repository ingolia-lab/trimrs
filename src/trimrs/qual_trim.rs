//! Quality trimming from either terminus of the read, start or end.
//!
//! The low-quality region can include some bases with quality above
//! the `threshold`. The algorithm computes a running sum of the
//! `qual[i] - threshold` difference starting from the terminus and
//! keeps track of the lowest score seen, along with its
//! position. This process stops as soon as the running sum is
//! positive, and the low-quality region then starts at the position
//! of this running sum minimum. When multiple minima are present, the
//! longest low-quality region is chosen.
//!
//! Since qualities above `threshold` contribute positive values to
//! this running sum, the low-quality region is guaranteed to start at
//! a low-quality (`<= threshold`) site preceded by a high-quality (`>
//! trheshold`) site, except in the trivial cases where the terminal
//! position is high-quality and no trimming is performed, or the
//! low-quality region encompasses the entire read.
//!
//! The algorithm is taken from `Cutadapt`, which in turn took it from
//! `BWA`.

/// Computes read length remaining after dropping low-quality bases
/// from the end.
///
/// The remaining length may range from `0` to `quals.len()`
/// inclusive.
///
/// # Arguments
/// * `threshold` is the quality threshold
/// * `quals` are the qualities
pub fn qual_cutoff_end(threshold: u8, quals: &[u8]) -> usize {
    let mut running_sum = 0;
    let mut lowest_score = 0;
    let mut lowest_offset = quals.len();
    let threshold_i = threshold as isize;

    for i in (0..quals.len()).rev() {
        running_sum += quals[i] as isize - threshold_i;

        if running_sum > 0 {
            break;
        } else if running_sum <= lowest_score {
            lowest_score = running_sum;
            lowest_offset = i;
        }
    }

    return lowest_offset;
}

pub fn qual_cutoff_start(threshold: u8, quals: &[u8]) -> usize {
    let mut running_sum = 0;
    let mut lowest_score = 0;
    let mut lowest_offset = 0;
    let threshold_i = threshold as isize;

    for i in 0..quals.len() {
        running_sum += quals[i] as isize - threshold_i;

        if running_sum > 0 {
            break;
        } else if running_sum <= lowest_score {
            lowest_score = running_sum;
            lowest_offset = i + 1;
        }
    }

    return lowest_offset;
}

// Cutadapt example for end trimming
// i  0     1   2   3    4    5    6    7    8   9
// q 42,   40, 26, 27,   8,   7,  11,   4,   2,  3
// Threshold 10
//   32,   30, 16, 17,  -2,  -3,   1,  -6,  -8, -7
// Running sums
//  (70), (38), 8, -8, -25, -23, -20, -21, -15, -7

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quality_trim_end() {
        assert_eq!(qual_cutoff_end(10, &vec![]), 0);
        
        assert_eq!(qual_cutoff_end(10, &vec![9]), 0);
        assert_eq!(qual_cutoff_end(10, &vec![10]), 0);
        assert_eq!(qual_cutoff_end(10, &vec![11]), 1);

        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 9]), 2);
        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 10]), 2);
        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 11]), 3);

        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 9, 20]), 4);
        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 10, 20]), 4);
        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 11, 20]), 4);

        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 10, 9, 10]), 2);
        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 11, 9, 10]), 3);
        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 11, 9, 9]), 3);

        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 9, 11, 9, 10]), 2);

        assert_eq!(qual_cutoff_end(10, &vec![20, 20, 9, 12, 9, 10]), 4);

        assert_eq!(qual_cutoff_end(10, &vec![42, 40, 26, 27, 8, 7, 11, 4, 2, 3]), 4);
    }

    #[test]
    fn quality_trim_start() {
        assert_eq!(qual_cutoff_start(10, &vec![]), 0);
        
        assert_eq!(qual_cutoff_start(10, &vec![9]), 1);
        assert_eq!(qual_cutoff_start(10, &vec![10]), 1);
        assert_eq!(qual_cutoff_start(10, &vec![11]), 0);

        assert_eq!(qual_cutoff_start(10, &vec![9, 20, 20]), 1);
        assert_eq!(qual_cutoff_start(10, &vec![10, 20, 20]), 1);
        assert_eq!(qual_cutoff_start(10, &vec![11, 20, 20]), 0);
        
        assert_eq!(qual_cutoff_start(10, &vec![20,  9, 20, 20]), 0);
        assert_eq!(qual_cutoff_start(10, &vec![20, 10, 20, 20]), 0);
        assert_eq!(qual_cutoff_start(10, &vec![20, 11, 20, 20]), 0);

        assert_eq!(qual_cutoff_start(10, &vec![10, 9, 10, 20, 20]), 3);
        assert_eq!(qual_cutoff_start(10, &vec![10, 9, 11, 20, 20]), 2);
        assert_eq!(qual_cutoff_start(10, &vec![ 9, 9, 11, 20, 20]), 2);

        assert_eq!(qual_cutoff_start(10, &vec![10, 9, 11, 9, 20, 20]), 4);

        assert_eq!(qual_cutoff_start(10, &vec![10, 9, 12, 9, 20, 20]), 2);

        
        assert_eq!(qual_cutoff_start(10, &vec![3, 2, 4, 11, 7, 8, 27, 26, 40, 42]), 6);
    }

}

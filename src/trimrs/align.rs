use std::default::Default;

use anyhow::{bail, ensure, Result};

use crate::cutadapt_encode::*;

#[derive(Clone,Copy,Debug,PartialEq,Eq,PartialOrd,Ord)]
enum Origin {
    // Origin positive
    QueryStart(usize),
    // Origin negative
    RefStart(usize),
}

impl std::default::Default for Origin {
    fn default() -> Self { Origin::RefStart(0) }
}

// Structure for a DP matrix entry
#[derive(Clone, Copy, Debug, Default)]
struct Entry {
    cost: usize,
    matches: usize,
    origin: Origin,
}

#[derive(Clone, Copy, Debug, Default)]
struct Match {
    origin: Origin,
    cost: usize,
    matches: usize,
    ref_stop: usize,
    query_stop: usize,
}

///     Representation of the dynamic-programming matrix.
///
///     This is used only when debugging is enabled in the Aligner class since the
///     matrix is normally not stored in full.
///
///     Entries in the matrix may be None, in which case that value was not
///     computed.
#[derive(Clone, Debug)]
pub struct DPMatrix {
    rows: Vec<Vec<Option<usize>>>,
    reference: Vec<u8>,
    query: Vec<u8>,
}

impl DPMatrix {
    pub fn new(reference: &[u8], query: &[u8]) -> Self {
        let m = reference.len();
        let n = query.len();
        let rows = vec![vec![None; n + 1]; m + 1];
        Self {
            rows: rows,
            reference: reference.to_vec(),
            query: query.to_vec(),
        }
    }

    /// Set an entry in the dynamic programming matrix.
    pub fn set_entry(&mut self, i: usize, j: usize, cost: usize) {
        self.rows[i][j] = Some(cost);
    }
}

impl std::fmt::Display for DPMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let headers = self
            .query
            .iter()
            .map(|c| format!("{:>2}", *c as char))
            .collect::<Vec<String>>();
        let header = "     ".to_string() + &headers.join(" ");
        let mut rows = vec![header];
        for (q, row) in std::iter::once(b' ')
            .chain(self.reference.iter().copied())
            .zip(self.rows.iter())
        {
            let cells = row
                .iter()
                .map(|v| match v {
                    None => "  ".to_string(),
                    Some(c) => format!("{:>2}", c),
                })
                .collect::<Vec<String>>();
            let r = format!("{} {}", q as char, cells.join(" "));

            rows.push(r);
        }
        write!(f, "{}", rows.join("\n"))
    }
}

#[derive(PartialEq,Eq,PartialOrd,Ord,Clone,Copy,Hash,Debug)]
pub enum AlignEnds {
    Global,
    LocalStart,
    LocalStop,
    Local
}

impl AlignEnds {
    pub fn start_local(self) -> bool {
        self == Self::LocalStart || self == Self::Local
    }

    pub fn stop_local(self) -> bool {
        self == Self::LocalStop || self == Self::Local
    }
}

#[derive(PartialEq,Eq,PartialOrd,Ord,Clone,Copy,Hash,Debug)]
pub enum AlignMatching {
    NoWildcard,
    RefWildcard,
    QueryWildcard
}

impl AlignMatching {
    pub fn ref_wildcard(self) -> bool {
        self == AlignMatching::RefWildcard
    }

    pub fn query_wildcard(self) -> bool {
        self == AlignMatching::QueryWildcard
    }
}

/// Parameters of an [`Aligner`](struct.Aligner.html) match. The
/// alignment itself is not computed, only the starting and stopping
/// positions, along with the number of matches and errors. Alignment
/// positions are reported in Rust coordinate conventions, with
/// _start_ as the first position in the alignment and _stop_ not
/// included, i.e., _start_ to (_stop_-1) inclusive.
#[derive(PartialEq,Eq,PartialOrd,Ord,Debug,Clone)]
pub struct Location {
    refstart: usize,
    refstop: usize,
    querystart: usize,
    querystop: usize,
    matches: usize,
    errors: usize,
}

impl Location {
    /// Starting position on reference sequence
    pub fn refstart(&self) -> usize {
        self.refstart
    }

    /// Stopping position on reference sequence, as a half-open
    /// interval, i.e., alignment does not include this position.
    pub fn refstop(&self) -> usize {
        self.refstop
    }

    /// Starting position on the query sequence
    pub fn querystart(&self) -> usize {
        self.querystart
    }

    /// Stopping position on the query sequence, as a half-open
    /// interval, i.e., alignment does not include this position.
    pub fn querystop(&self) -> usize {
        self.querystop
    }

    pub fn matches(&self) -> usize {
        self.matches
    }

    pub fn errors(&self) -> usize {
        self.errors
    }
}

/// Configuration structure for `Aligner`.
///
/// Alignment parameters are named fields in the structure.
#[derive(Clone,Debug,PartialEq,PartialOrd)]
pub struct AlignerConf {
    /// Maximum error rate
    pub max_error_rate: f64,

    /// Start and/or end gaps in the reference sequence
    pub reference_ends: AlignEnds,

    /// Start and/or end gaps in the query sequence
    pub query_ends: AlignEnds,

    /// Use IUPAC wildcards in reference or query
    pub matching: AlignMatching,

    /// Cost of insertion or deletion
    pub indel_cost: usize,

    /// Minimum overlap to report a match
    pub min_overlap: usize,
}
    

///    Find a full or partial occurrence of a query string in a reference string
///    allowing errors (mismatches, insertions, deletions).
///
///
///    Mismatches are counted as one error, and insertions and
///    deletions are counted as one or mor errors.
///
///    Semi-global alignments allow skipping a suffix and/or prefix of query or
///    reference at no cost. Combining semi-global alignment with edit distance is
///    a bit unusual because the trivial “optimal” solution at edit distance 0
///    would be to skip all of the reference and all of the query, like this:
///
///    ```
///        REFERENCE-----
///        ---------QUERY
///    ```
///
///    Conceptually, the algorithm used here instead tests all possible overlaps
///    between the two sequences and chooses the overlap which maximizes the
///    number of matches in the overlapping part, while the error rate must not
///    go above a threshold.
///
///    The `reference_ends` and `query_ends` parameters control
///    whether prefixes and suffixes are skipped. If both are set to
///    `LocalStart` or `Local`, allowing prefixes to be skipped at no
///    cost, then a prefix of _either_ `reference` _or_ `query` is
///    skipped, never both. Similarly, if both reference and query
///    allow local alignment at the stop, then only one string will
///    have a skipped suffix.
///
///    The skipped parts are described with two intervals (_refstart_,
///    _refstop_) and (_querystart_, _querystop_).
///
///    For example, an optimal semiglobal alignment of SISSI and MISSISSIPPI looks like this:
///
///    ```
///    ref   ---SISSI---
///    query MISSISSIPPI
///    ```
///
///    ```
///    refstart, refstop = 0, 5
///    querystart, querystop = 3, 8
///    (with zero errors)
///    ```
///
///    The aligned parts are `reference[refstart..refstop]` and
///    `query[querystart..querystop]`.
///
///    The error rate is: _errors_ / _length_ where _length_ is
///    `(refstop - refstart)`.
///
///    An optimal alignment fulfills all of these criteria:
///
///    - its error rate is at most `max_error_rate`
///    - Among those alignments with an error rate no greater than `max_error_rate`, the alignment contains
///      a maximal number of matches (there is no alignment with more matches).
///    - If there are multiple alignments with the same number of matches, then one that
///      has minimal number of errors is chosen.
///    - If there are still multiple candidates, choose the alignment that starts at the
///      leftmost position within the read.
///
///    The alignment itself is not returned, only the information in
///    [`Location`](struct.Location.html).
///    It is always the case that at least one of `refstart` and `querystart` is zero.
///
///    IUPAC wildcard characters can be allowed in the reference or
///    the query by setting the appropriate flags. All non-IUPAC
///    characters compare as non-equal, and when IUPAC wildcards are
///    not enabled, anything except `A`, `C`, `G`, `T`, and `U`
///    compares as non-equal to everything.
#[derive(Clone, Debug)]
pub struct Aligner {
    column: Vec<Entry>,
    max_error_rate: f64,
    reference_ends: AlignEnds,
    query_ends: AlignEnds,
    insertion_cost: usize,
    deletion_cost: usize,
    min_overlap: usize,
    matching: AlignMatching,
    debug: bool,
    dpmatrix: Option<DPMatrix>,
    reference: Vec<u8>,
    breference: Vec<u8>,
    effective_length: isize,
    n_counts: Vec<isize>,
    bquery: Vec<u8>,
}

const INIT_QUERY_LEN: usize = 256;

impl Aligner {
    /// Creates a new aligner with a specified alignment configuration
    /// and reference sequence.
    pub fn new(
        conf: &AlignerConf,
        reference: &[u8],
    ) -> Result<Self> {
        let m = reference.len();
        let mut n_counts = vec![0; m + 1];
        let mut effective_length = m as isize;
        let mut breference = reference.to_vec();

        let mut n_count = 0;
        for i in 0..m {
            n_counts[i] = n_count;
            if reference[i] == b'n' || reference[i] == b'N' {
                n_count += 1;
            }
        }
        *n_counts.last_mut().unwrap() = n_count;
        assert_eq!(
            n_counts[m],
            reference
                .iter()
                .copied()
                .filter(|&c| c == b'N' || c == b'n')
                .count() as isize
        );
        match conf.matching {
            AlignMatching::RefWildcard => {
                effective_length = m as isize - n_counts[m];
                if effective_length == 0 {
                    bail!("Cannot have only N wildcards in the sequence");
                }
                encode_iupac_vec(reference, &mut breference);
            },
            AlignMatching::NoWildcard |
            AlignMatching::QueryWildcard => {
                encode_acgt_vec(reference, &mut breference);
            },
        };

        ensure!(conf.indel_cost >= 1, "Indel cost must be at least 1");
        ensure!(conf.min_overlap >= 1, "Min overlap must be at least 1");

        Ok(Aligner {
            column: vec![Entry::default(); m + 1],
            max_error_rate: conf.max_error_rate,
            reference_ends: conf.reference_ends,
            query_ends: conf.query_ends,
            matching: conf.matching,
            min_overlap: conf.min_overlap,
            debug: false,
            dpmatrix: None,
            reference: reference.to_vec(),
            breference: breference,
            effective_length: effective_length,
            n_counts: n_counts,
            insertion_cost: conf.indel_cost,
            deletion_cost: conf.indel_cost,
            bquery: Vec::with_capacity(INIT_QUERY_LEN),
        })
    }

    /// Returns the length of the reference sequence
    pub fn m(&self) -> usize {
        self.reference.len()
    }

    /// Returns the dynamic programming matrix, which is `None` unless
    /// debugging has been enabled.
    pub fn dpmatrix(&self) -> &Option<DPMatrix> {
        &self.dpmatrix
    }

    /// Store the dynamic programming matrix while running
    /// [`locate()`](#method.locate) and make it available by
    /// [`dpmatrix()`](#method.dpmatrix).
    pub fn enable_debug(&mut self) {
        self.debug = true;
    }

    //         locate(query) -> (refstart, refstop, querystart, querystop, matches, errors)

    //         Find the query within the reference associated with this aligner. The
    //         intervals (querystart, querystop) and (refstart, refstop) give the
    //         location of the match.

    //         That is, the substrings query[querystart:querystop] and
    //         self.reference[refstart:refstop] were found to align best to each other,
    //         with the given number of matches and the given number of errors.

    //         The alignment itself is not returned.
    pub fn locate(&mut self, query: &[u8]) -> Option<Location> {
        let s1 = &self.breference;
        let m = self.m();
        let n = query.len();
        let column = &mut self.column;
        let max_error_rate = self.max_error_rate;
        let stop_in_query = self.query_ends.stop_local();

        match self.matching {
            AlignMatching::QueryWildcard => {
                encode_iupac_vec(query, &mut self.bquery);
            },
            AlignMatching::NoWildcard |
            AlignMatching::RefWildcard => {
                encode_acgt_vec(query, &mut self.bquery);
            },
        };
        let s2 = &self.bquery;
        //         DP Matrix:
        //                    query (j)
        //                  ----------> n
        //                 |
        //         ref (i) |
        //                 |
        //                 V
        //                m

        // # maximum no. of errors
        let k = (max_error_rate * m as f64).floor() as usize;

        // # Determine largest and smallest column we need to compute
        let max_n = if !self.query_ends.start_local() {
            // # costs can only get worse after column m
            usize::min(n, m + k)
        } else {
            n
        };
        let min_n = if !self.query_ends.stop_local() {
            if n > m + k {
                n - m - k
            } else {
                0
            }
        } else {
            0
        };

        // # Fill column min_n.
        // #
        // # Four cases:
        // # not startin1, not startin2: c(i,j) = max(i,j); origin(i, j) = 0
        // #     startin1, not startin2: c(i,j) = j       ; origin(i, j) = min(0, j - i)
        // # not startin1,     startin2: c(i,j) = i       ; origin(i, j) =
        // #     startin1,     startin2: c(i,j) = min(i,j)

        // # TODO (later)
        // # fill out columns only until 'last'
        if !self.reference_ends.start_local() && !self.query_ends.start_local() {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = usize::max(i, min_n) * self.insertion_cost;
                column[i].origin = Origin::RefStart(0);
            }
        } else if self.reference_ends.start_local() && !self.query_ends.start_local() {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = min_n * self.insertion_cost;
                column[i].origin = Origin::RefStart(if i > min_n {
                    i - min_n
                } else {
                    0
                });
                    //                column[i].origin = isize::min(0, min_n - i as isize);
                    // max(0, i - min_n)
            }
        } else if !self.reference_ends.start_local() && self.query_ends.start_local() {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = i * self.insertion_cost;
                column[i].origin = Origin::QueryStart(if min_n > i {
                    min_n - i
                } else {
                    0
                });
            }
        } else {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = usize::min(i, min_n) * self.insertion_cost;
                column[i].origin = if min_n > i {
                    Origin::QueryStart(min_n - i)
                } else {
                    Origin::RefStart(i - min_n)
                };
            }
        }

        if self.debug {
            let mut dpmatrix = DPMatrix::new(&self.reference, query);
            for i in 0..(m + 1) {
                dpmatrix.set_entry(i, min_n, column[i].cost);
            }
            self.dpmatrix = Some(dpmatrix);
        }

        let mut best = Match::default();
        best.ref_stop = m;
        best.query_stop = n;
        best.cost = m + n;
        best.origin = Origin::RefStart(0);
        best.matches = 0;

        // # Ukkonen's trick: index of the last cell that is at most k
        let mut last = isize::min(m as isize, k as isize + 1);
        if self.reference_ends.start_local() {
            last = m as isize;
        }

        // # We keep only a single column of the DP matrix in memory.
        // # To access the diagonal cell to the upper left,
        // # we store it here before overwriting it.
        let mut diag_entry;

        // # iterate over columns
        for j in (min_n + 1)..(max_n + 1) {
            // # remember first entry before overwriting
            diag_entry = column[0];

            // # fill in first entry in this column
            if self.query_ends.start_local() {
                column[0].origin = Origin::QueryStart(j);
            } else {
                column[0].cost = j * self.insertion_cost;
            }

            for i in 0..(last as usize) {
                let origin;
                let cost;
                let matches;

                let characters_equal = (s1[i] & s2[(j - 1)]) != 0;

                if characters_equal {
                    // # If the characters match, skip computing costs for
                    // # insertion and deletion as they are at least as high.
                    cost = diag_entry.cost;
                    origin = diag_entry.origin;
                    matches = diag_entry.matches + 1;
                } else {
                    // # Characters do not match.
                    let cost_diag = diag_entry.cost + 1;
                    let cost_deletion = column[i + 1].cost + self.deletion_cost;
                    let cost_insertion = column[i].cost + self.insertion_cost;

                    if cost_diag <= cost_deletion && cost_diag <= cost_insertion {
                        // # MISMATCH
                        cost = cost_diag;
                        origin = diag_entry.origin;
                        matches = diag_entry.matches;
                    } else if cost_insertion < cost_deletion {
                        // # INSERTION
                        cost = cost_insertion;
                        origin = column[i].origin;
                        matches = column[i].matches;
                    } else {
                        // # DELETION
                        cost = cost_deletion;
                        origin = column[i + 1].origin;
                        matches = column[i + 1].matches;
                    }
                }

                // # Remember the current cell for next iteration
                diag_entry = column[i + 1];

                column[i + 1].cost = cost;
                column[i + 1].origin = origin;
                column[i + 1].matches = matches;
            }
            if let Some(dpmatrix) = self.dpmatrix.as_mut() {
                for i in 0..(last as usize + 1) {
                    dpmatrix.set_entry(i, j, column[i].cost);
                }
            }
            while last >= 0 && column[last as usize].cost > k {
                last -= 1;
            }

            // # last can be -1 here, but will be incremented next.
            // # TODO if last is -1, can we stop searching?
            if last < m as isize {
               last += 1;
            } else if self.query_ends.stop_local() {
                // # Found a match. If requested, find best match in last row.
                // # length of the aligned part of the reference
                //                let length = m as isize + isize::min(column[m].origin, 0);
                let length = match column[m].origin {
                    Origin::QueryStart(_) => m as isize,
                    Origin::RefStart(r) => (m - r) as isize,
                };
                let cur_effective_length = if self.matching.ref_wildcard() {
                    if length < m as isize {
                        // # Recompute effective length so that it only takes into
                        // # account the matching suffix of the reference
                        length - self.n_counts[length as usize]
                    } else {
                        self.effective_length
                    }
                } else {
                    length
                };
                let cost = column[m].cost;
                let matches = column[m].matches;
                if length >= self.min_overlap as isize
                    && (cost as f64) <= (cur_effective_length as f64) * max_error_rate
                    && (matches > best.matches || (matches == best.matches && cost < best.cost))
                {
                    // # update
                    best.matches = matches;
                    best.cost = cost;
                    best.origin = column[m].origin;
                    best.ref_stop = m;
                    best.query_stop = j;
                    if cost == 0 && matches == m {
                        // # exact match, stop early
                        break;
                    }
                    // # column finished
                }
            }
        }

        if max_n == n {
            let first_i = if self.reference_ends.stop_local() { 0 } else { m };

            // # search in last column # TODO last?
            for i in first_i..(m + 1) {
                let length = match column[i].origin {
                    Origin::QueryStart(_) => i as isize,
                    Origin::RefStart(s) => (i - s) as isize,
                };
                let cost = column[i].cost;
                let matches = column[i].matches;
                let cur_effective_length = if self.matching.ref_wildcard() {
                    if length < m as isize {
                        // # Recompute effective length so that it only takes into
                        // # account the matching part of the reference
                        let ref_start = match column[i].origin {
                            Origin::QueryStart(_) => 0,
                            Origin::RefStart(r) => r,
                        };
                        assert!(0 <= ref_start);
                        assert!(ref_start <= m);
                        length - (self.n_counts[i] - self.n_counts[ref_start])
                    } else {
                        self.effective_length
                    }
                } else {
                    length
                };

                assert!(0 <= cur_effective_length);
                assert!(cur_effective_length <= length);
                assert!(cur_effective_length <= self.effective_length);

                if length >= self.min_overlap as isize
                    && (cost as f64) <= cur_effective_length as f64 * max_error_rate
                    && (matches > best.matches || (matches == best.matches && cost < best.cost))
                {
                    // # update best
                    best.matches = matches;
                    best.cost = cost;
                    best.origin = column[i].origin;
                    best.ref_stop = i;
                    best.query_stop = n;
                }
            }
        }

        if best.cost == m + n {
            // # best.cost was initialized with this value.
            // # If it is unchanged, no alignment was found that has
            // # an error rate within the allowed range.
            return None;
        }

        let refstart;
        let querystart;
        match best.origin {
            Origin::QueryStart(start) => {
                refstart = 0;
                querystart = start;
            },
            Origin::RefStart(start) => {
                refstart = start;
                querystart = 0;
            },
        };
        assert!(best.ref_stop > refstart); // # Do not return empty alignments.
        return Some(
            Location {
                refstart: refstart,
                refstop: best.ref_stop,
                querystart: querystart,
                querystop: best.query_stop,
                matches: best.matches,
                errors: best.cost,
            });
    }
}

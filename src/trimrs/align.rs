use std::default::Default;

use anyhow::{bail, ensure, Result};

use crate::cutadapt_encode::*;

// Structure for a DP matrix entry
#[derive(Clone, Copy, Debug, Default)]
struct Entry {
    cost: isize,
    matches: isize,
    origin: isize,
}

#[derive(Clone, Copy, Debug, Default)]
struct Match {
    origin: isize,
    cost: isize,
    matches: isize,
    ref_stop: isize,
    query_stop: isize,
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
    rows: Vec<Vec<Option<isize>>>,
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
    pub fn set_entry(&mut self, i: isize, j: isize, cost: isize) {
        self.rows[i as usize][j as usize] = Some(cost);
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

///    Find a full or partial occurrence of a query string in a reference string
///    allowing errors (mismatches, insertions, deletions).

///    By default, unit costs are used, meaning that mismatches, insertions and
///    deletions are counted as one error (edit distance).

///    Semi-global alignments allow skipping a suffix and/or prefix of query or
///    reference at no cost. Combining semi-global alignment with edit distance is
///    a bit unusual because the trivial “optimal” solution at edit distance 0
///    would be to skip all of the reference and all of the query, like this:

///        REFERENCE-----
///        ---------QUERY

///    Conceptually, the algorithm used here instead tests all possible overlaps
///    between the two sequences and chooses the overlap which maximizes the
///    number of matches in the overlapping part, while the error rate must not
///    go above a threshold.

///    TODO working here

///    To allow skipping of a prefix of string1 at no cost, set the
///    START_IN_REFERENCE flag.
///    To allow skipping of a prefix of string2 at no cost, set the
///    START_IN_QUERY flag.
///    If both are set, a prefix of string1 or of string1 is skipped,
///    never both.
///    Similarly, set STOP_IN_REFERENCE and STOP_IN_QUERY to
///    allow skipping of suffixes of string1 or string2. Again, when both
///    flags are set, never suffixes in both strings are skipped.
///    If all flags are set, this results in standard semiglobal alignment.

///    The skipped parts are described with two intervals (start1, stop1),
///    (start2, stop2).

///    For example, an optimal semiglobal alignment of SISSI and MISSISSIPPI looks like this:

///    ```
///    ---SISSI---
///    MISSISSIPPI
///    ```

///    start1, stop1 = 0, 5
///    start2, stop2 = 3, 8
///    (with zero errors)

///    The aligned parts are string1[start1:stop1] and string2[start2:stop2].

///    The error rate is: errors / length where length is (stop1 - start1).

///    An optimal alignment fulfills all of these criteria:

///    - its error_rate is at most max_error_rate
///    - Among those alignments with error_rate <= max_error_rate, the alignment contains
///      a maximal number of matches (there is no alignment with more matches).
///    - If there are multiple alignments with the same no. of matches, then one that
///      has minimal no. of errors is chosen.
///    - If there are still multiple candidates, choose the alignment that starts at the
///      leftmost position within the read.

///    The alignment itself is not returned, only the tuple
///    (start1, stop1, start2, stop2, matches, errors), where the first four fields have the
///    meaning as described, matches is the number of matches and errors is the number of
///    errors in the alignment.

///    It is always the case that at least one of start1 and start2 is zero.

///    IUPAC wildcard characters can be allowed in the reference and the query
///    by setting the appropriate flags.

///    If neither flag is set, the full ASCII alphabet is used for comparison.
///    If any of the flags is set, all non-IUPAC characters in the sequences
///    compare as 'not equal'.
#[derive(Clone, Debug)]
pub struct Aligner {
    column: Vec<Entry>,
    max_error_rate: f64,
    start_in_reference: bool,
    start_in_query: bool,
    stop_in_reference: bool,
    stop_in_query: bool,
    insertion_cost: isize,
    deletion_cost: isize,
    min_overlap: isize,
    wildcard_ref: bool,
    wildcard_query: bool,
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
    pub fn m(&self) -> usize {
        self.reference.len()
    }

    pub fn new(
        reference: &[u8],
        max_error_rate: f64,
        flags: isize,
        wildcard_ref: bool,
        wildcard_query: bool,
        indel_cost: isize,
        min_overlap: isize,
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
            *n_counts.last_mut().unwrap() = n_count;
            assert!(
                n_counts[m]
                    == reference
                        .iter()
                        .copied()
                        .filter(|&c| c == b'N' || c == b'n')
                        .count() as isize
            );
            if wildcard_ref {
                effective_length = m as isize - n_counts[m];
                if effective_length == 0 {
                    bail!("Cannot have only N wildcards in the sequence");
                }
                encode_iupac_vec(reference, &mut breference);
            } else if wildcard_query {
                encode_acgt_vec(reference, &mut breference);
            } else {
                breference.make_ascii_uppercase();
            }
        }

        ensure!(indel_cost >= 1, "Indel cost must be at least 1");
        ensure!(min_overlap >= 1, "Min overlap must be at least 1");

        Ok(Aligner {
            column: vec![Entry::default(); m + 1],
            max_error_rate: max_error_rate,
            start_in_reference: flags & 1 > 0,
            start_in_query: flags & 2 > 0,
            stop_in_reference: flags & 4 > 0,
            stop_in_query: flags & 8 > 0,
            wildcard_ref: wildcard_ref,
            wildcard_query: wildcard_query,
            min_overlap: min_overlap,
            debug: false,
            dpmatrix: None,
            reference: reference.to_vec(),
            breference: breference,
            effective_length: effective_length,
            n_counts: n_counts,
            insertion_cost: indel_cost,
            deletion_cost: indel_cost,
            bquery: Vec::with_capacity(INIT_QUERY_LEN),
        })
    }

    //     property dpmatrix:
    //         The dynamic programming matrix as a DPMatrix object. This attribute is
    //         usually None, unless debugging has been enabled with enable_debug().
    pub fn dpmatrix(&self) -> &Option<DPMatrix> {
        &self.dpmatrix
    }

    ///         Store the dynamic programming matrix while running the locate() method
    ///         and make it available in the .dpmatrix attribute.
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
    pub fn locate(&mut self, query: &[u8]) -> Option<(isize, isize, isize, isize, isize, isize)> {
        let s1 = &self.reference;
        let m = self.m();
        let n = query.len() as isize;
        let column = &mut self.column;
        let max_error_rate = self.max_error_rate;
        let stop_in_query = self.stop_in_query;

        let compare_ascii = if self.wildcard_query {
            encode_iupac_vec(query, &mut self.bquery);
            false
        } else if self.wildcard_ref {
            encode_acgt_vec(query, &mut self.bquery);
            false
        } else {
            self.bquery.clear();
            self.bquery.extend(query.iter());
            self.bquery.make_ascii_uppercase();
            true
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
        let k = (max_error_rate * m as f64).floor() as isize;

        // # Determine largest and smallest column we need to compute
        let max_n = if !self.start_in_query {
            // # costs can only get worse after column m
            isize::min(n, m as isize + k)
        } else {
            n
        };
        let min_n = if !self.stop_in_query {
            isize::max(0, n - m as isize - k)
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
        if !self.start_in_reference && !self.start_in_query {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = isize::max(i as isize, min_n) * self.insertion_cost;
                column[i].origin = 0;
            }
        } else if self.start_in_reference && !self.start_in_query {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = min_n * self.insertion_cost;
                column[i].origin = isize::min(0, min_n - i as isize);
            }
        } else if !self.start_in_reference && self.start_in_query {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = i as isize * self.insertion_cost;
                column[i].origin = isize::max(0, min_n - i as isize);
            }
        } else {
            for i in 0..(m + 1) {
                column[i].matches = 0;
                column[i].cost = isize::min(i as isize, min_n) * self.insertion_cost;
                column[i].origin = min_n - i as isize;
            }
        }

        if self.debug {
            let mut dpmatrix = DPMatrix::new(&self.reference, query);
            for i in 0..(m + 1) {
                dpmatrix.set_entry(i as isize, min_n, column[i as usize].cost);
            }
            self.dpmatrix = Some(dpmatrix);
        }

        let mut best = Match::default();
        best.ref_stop = m as isize;
        best.query_stop = n;
        best.cost = m as isize + n;
        best.origin = 0;
        best.matches = 0;

        // # Ukkonen's trick: index of the last cell that is at most k
        let mut last = isize::min(m as isize, k + 1);
        if self.start_in_reference {
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
            if self.start_in_query {
                column[0].origin = j;
            } else {
                column[0].cost = j * self.insertion_cost;
            }

            for i in 1..(last + 1) {
                let origin;
                let cost;
                let matches;

                let characters_equal = if compare_ascii {
                    (s1[(i - 1) as usize] == s2[(j - 1) as usize])
                } else {
                    (s1[(i - 1) as usize] & s2[(j - 1) as usize]) != 0
                };

                if characters_equal {
                    // # If the characters match, skip computing costs for
                    // # insertion and deletion as they are at least as high.
                    cost = diag_entry.cost;
                    origin = diag_entry.origin;
                    matches = diag_entry.matches + 1;
                } else {
                    // # Characters do not match.
                    let cost_diag = diag_entry.cost + 1;
                    let cost_deletion = column[i as usize].cost + self.deletion_cost;
                    let cost_insertion = column[(i - 1) as usize].cost + self.insertion_cost;

                    if cost_diag <= cost_deletion && cost_diag <= cost_insertion {
                        // # MISMATCH
                        cost = cost_diag;
                        origin = diag_entry.origin;
                        matches = diag_entry.matches;
                    } else if cost_insertion < cost_deletion {
                        // # INSERTION
                        cost = cost_insertion;
                        origin = column[(i - 1) as usize].origin;
                        matches = column[(i - 1) as usize].matches;
                    } else {
                        // # DELETION
                        cost = cost_deletion;
                        origin = column[i as usize].origin;
                        matches = column[i as usize].matches;
                    }
                }

                // # Remember the current cell for next iteration
                diag_entry = column[i as usize];

                column[i as usize].cost = cost;
                column[i as usize].origin = origin;
                column[i as usize].matches = matches;
            }
            if let Some(dpmatrix) = self.dpmatrix.as_mut() {
                for i in 0..(last + 1) {
                    dpmatrix.set_entry(i, j, column[i as usize].cost);
                }
            }
            while last >= 0 && column[last as usize].cost > k {
                last -= 1;
            }

            // # last can be -1 here, but will be incremented next.
            // # TODO if last is -1, can we stop searching?
            if last < m as isize {
                last += 1;
            } else if stop_in_query {
                // # Found a match. If requested, find best match in last row.
                // # length of the aligned part of the reference
                let length = m as isize + isize::min(column[m].origin, 0);
                let cur_effective_length = if self.wildcard_ref {
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
                let cost = column[m as usize].cost;
                let matches = column[m as usize].matches;
                if length >= self.min_overlap
                    && (cost as f64) <= (cur_effective_length as f64) * max_error_rate
                    && (matches > best.matches || (matches == best.matches && cost < best.cost))
                {
                    // # update
                    best.matches = matches;
                    best.cost = cost;
                    best.origin = column[m as usize].origin;
                    best.ref_stop = m as isize;
                    best.query_stop = j;
                    if cost == 0 && matches == m as isize {
                        // # exact match, stop early
                        break;
                    }
                    // # column finished
                }
            }
        }

        if max_n == n {
            let first_i = if self.stop_in_reference { 0 } else { m };

            // # search in last column # TODO last?
            for i in first_i..(m + 1) {
                let length = i as isize + isize::min(column[i].origin, 0);
                let cost = column[i].cost;
                let matches = column[i].matches;
                let cur_effective_length = if self.wildcard_ref {
                    if length < m as isize {
                        // # Recompute effective length so that it only takes into
                        // # account the matching part of the reference
                        let ref_start = -isize::min(column[i].origin, 0);
                        assert!(0 <= ref_start);
                        assert!(ref_start <= m as isize);
                        length - (self.n_counts[i] - self.n_counts[ref_start as usize])
                    } else {
                        self.effective_length
                    }
                } else {
                    length
                };

                assert!(0 <= cur_effective_length);
                assert!(cur_effective_length <= length);
                assert!(cur_effective_length <= self.effective_length);

                if length >= self.min_overlap
                    && (cost as f64) <= cur_effective_length as f64 * max_error_rate
                    && (matches > best.matches || (matches == best.matches && cost < best.cost))
                {
                    // # update best
                    best.matches = matches;
                    best.cost = cost;
                    best.origin = column[i].origin;
                    best.ref_stop = i as isize;
                    best.query_stop = n;
                }
            }
        }

        if best.cost == m as isize + n {
            // # best.cost was initialized with this value.
            // # If it is unchanged, no alignment was found that has
            // # an error rate within the allowed range.
            return None;
        }

        let start1;
        let start2;
        if best.origin >= 0 {
            start1 = 0;
            start2 = best.origin;
        } else {
            start1 = -best.origin;
            start2 = 0;
        }
        assert!(best.ref_stop - start1 > 0); // # Do not return empty alignments.
        return Some((
            start1,
            best.ref_stop,
            start2,
            best.query_stop,
            best.matches,
            best.cost,
        ));
    }
}

//     def __dealloc__(self):
//         PyMem_Free(self.column)
//         PyMem_Free(self.n_counts)

// cdef class PrefixComparer:
//     """
//     A version of the Aligner that is specialized in the following way:

//     - it does not allow indels
//     - it allows only 5' anchored adapters

//     This is a separate class, not simply a function, in order to be able
//     to cache the reference (avoiding to convert it from str to bytes on
//     every invocation)
//     """
//     cdef:
//         bytes reference
//         bint wildcard_ref
//         bint wildcard_query
//         int m
//         int max_k  # max. number of errors
//         readonly int effective_length
//         int min_overlap

//     # __init__ instead of __cinit__ because we need to override this in SuffixComparer
//     def __init__(
//         self,
//         str reference,
//         double max_error_rate,
//         bint wildcard_ref=False,
//         bint wildcard_query=False,
//         int min_overlap=1,
//     ):
//         self.wildcard_ref = wildcard_ref
//         self.wildcard_query = wildcard_query
//         self.m = len(reference)
//         self.effective_length = self.m
//         if self.wildcard_ref:
//             self.effective_length -= reference.count('N') - reference.count('n')
//             if self.effective_length == 0:
//                 raise ValueError("Cannot have only N wildcards in the sequence")
//         if not (0 <= max_error_rate <= 1.):
//             raise ValueError("max_error_rate must be between 0 and 1")
//         self.max_k = int(max_error_rate * self.effective_length)
//         self.reference = reference.encode('ascii').upper()
//         if min_overlap < 1:
//             raise ValueError("min_overlap must be at least 1")
//         self.min_overlap = min_overlap
//         if self.wildcard_ref:
//             self.reference = self.reference.translate(IUPAC_TABLE)
//         elif self.wildcard_query:
//             self.reference = self.reference.translate(ACGT_TABLE)

//     def __repr__(self):
//         return "{}(reference={!r}, max_k={}, wildcard_ref={}, "\
//             "wildcard_query={})".format(
//                 self.__class__.__name__,
//                 self.reference, self.max_k, self.wildcard_ref,
//                 self.wildcard_query)

//     def locate(self, str query):
//         """
//         Find out whether one string is the prefix of the other one, allowing
//         IUPAC wildcards in ref and/or query if the appropriate flag is set.

//         This is used to find an anchored 5' adapter (type 'FRONT') in the 'no indels' mode.
//         This is very simple as only the number of errors needs to be counted.

//         This function returns a tuple compatible with what Aligner.locate outputs.
//         """
//         cdef:
//             bytes query_bytes = query.encode('ascii')
//             char* r_ptr = self.reference
//             char* q_ptr
//             int i, matches = 0
//             int n = len(query_bytes)
//             int length = min(self.m, n)
//             bint compare_ascii = False
//             int errors

//         if self.wildcard_query:
//             query_bytes = query_bytes.translate(IUPAC_TABLE)
//         elif self.wildcard_ref:
//             query_bytes = query_bytes.translate(ACGT_TABLE)
//         else:
//             query_bytes = query_bytes.upper()
//             compare_ascii = True
//         q_ptr = query_bytes

//         if compare_ascii:
//             for i in range(length):
//                 if r_ptr[i] == q_ptr[i]:
//                     matches += 1
//         else:
//             for i in range(length):
//                 if (r_ptr[i] & q_ptr[i]) != 0:
//                     matches += 1

//         errors = length - matches
//         if errors > self.max_k or length < self.min_overlap:
//             return None
//         return (0, length, 0, length, matches, length - matches)

// cdef class SuffixComparer(PrefixComparer):

//     def __init__(
//         self,
//         str reference,
//         double max_error_rate,
//         bint wildcard_ref=False,
//         bint wildcard_query=False,
//         int min_overlap=1,
//     ):
//         super().__init__(reference[::-1], max_error_rate, wildcard_ref, wildcard_query, min_overlap)

//     def locate(self, str query):
//         cdef int n = len(query)
//         result = super().locate(query[::-1])
//         if result is None:
//             return None
//         _, length, _, _, matches, errors = result
//         return (self.m - length, self.m, n - length, n, matches, errors)

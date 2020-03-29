use std::default::Default;

use trimrs::align;
use trimrs::encode::{encode_acgt_vec, encode_iupac_vec};

#[rustfmt::skip]
#[test]
fn basic_local_align() {
    test_align_both(b"CACCA", b"GACCACCATTA",
                    0.00, true, true, true, true, false, false, 1,
                    Some((0, 5, 3, 8, 5, 0)));
}

#[test]
fn bounds_align() {
    // start_in_query = false
    test_align_both(b"CACCA", b"GACCACCATTA",
                    0.00, true, false, true, true, false, false, 1,
                    None);
    test_align_both(b"GACCA", b"GACCACCATTA",
                    0.00, true, false, true, true, false, false, 1,
                    Some((0,5,0,5,5,0)));
    test_align_both(b"CACCA", b"GACCACCATTA",
                    0.25, true, false, true, true, false, false, 1,
                    Some((0,5,0,5,4,1)));

    // stop_in_query = false
    test_align_both(b"CACCA", b"GACCACCAATTA",
                    0.00, true, true, true, false, false, false, 1,
                    None);
    test_align_both(b"CACCA", b"ATTAGACCACCA",
                    0.00, true, true, true, false, false, false, 1,
                    Some((0,5,7,12,5,0)));
    test_align_both(b"GACCA", b"ATTAGACCACCA",
                    0.00, true, true, true, false, false, false, 1,
                    None);
    test_align_both(b"GACCA", b"ATTAGACCACCA",
                    0.25, true, true, true, false, false, false, 1,
                    Some((0,5,7,12,4,1)));

    // start_in_query and stop_in_query = false
    test_align_both(b"CACCA", b"GCACCAT", 0.00, true, false, true, false, false, false, 1, None);
    test_align_both(b"CACCA", b"GCACCA", 0.00, true, false, true, false, false, false, 1, None);
    test_align_both(b"CACCA", b"CACCAT", 0.00, true, false, true, false, false, false, 1, None);
    test_align_both(b"CACCA", b"GCACCAT", 0.25, true, false, true, false, false, false, 1, None);
    test_align_both(b"CACCA", b"GCACCA", 0.25, true, false, true, false, false, false, 1,
                    Some((0, 5, 0, 6, 5, 1)));
    test_align_both(b"CACCA", b"CACCAT", 0.25, true, false, true, false, false, false, 1,
                    Some((0, 5, 0, 6, 5, 1)));
}

#[test]
fn mismatch_align() {
    // Mismatches
    test_align_both(b"CACCA", b"GGACCAT", 0.00, true, true, true, true, false, false, 1, None);
    test_align_both(b"CACCA", b"GGACCAT", 0.25, true, true, true, true, false, false, 1,
                    Some((0, 5, 1, 6, 4, 1)));
    test_align_both(b"CACCA", b"GCTCCAT", 0.25, true, true, true, true, false, false, 1,
                    Some((0, 5, 1, 6, 4, 1)));
    test_align_both(b"CACCA", b"GCTCCAT", 0.25, true, true, true, true, false, false, 1,
                    Some((0, 5, 1, 6, 4, 1)));
    test_align_both(b"CACCA", b"GCAGCAT", 0.25, true, true, true, true, false, false, 1,
                    Some((0, 5, 1, 6, 4, 1)));
    test_align_both(b"CACCA", b"GCACGAT", 0.25, true, true, true, true, false, false, 1,
                    Some((0, 5, 1, 6, 4, 1)));
    // Apparently not biased towards match over indel at right-hand edge
    test_align_both(b"CACCA", b"GCACCTT", 0.25, true, true, true, true, false, false, 1,
                    Some((0, 5, 1, 5, 4, 1)));
}

#[test]
fn wildcards_align() {
    // Wildcards
    test_align_both(b"CACNA", b"GTGCACCATGT", 0.00, true, true, true, true, false, false, 1, None);
    test_align_both(b"CACNA", b"GTGCACCATGT", 0.25, true, true, true, true, false, false, 1,
                    Some((0, 5, 3, 8, 4, 1)));
    test_align_both(b"CACNA", b"GTGCACCATGT", 0.00, true, true, true, true, true, false, 1,
                    Some((0, 5, 3, 8, 5, 0)));    
    test_align_both(b"CACNA", b"GTGCACCATGT", 0.00, true, true, true, true, false, true, 1, None);

    // IUPAC degeneracy
    //   Reference A
    test_align_both(b"CACACAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACBCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACCCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACDCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACGCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACHCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACTCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACUCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACVCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACRCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACYCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACWCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACSCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACMCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACKCAC", b"GTGCACACACTGT", 0.0, true, true, true, true, true, false, 1, None);

    //   Reference C
    test_align_both(b"CACBCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACACAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACDCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACCCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACGCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACHCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACTCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACUCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACVCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACYCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACRCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACSCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACWCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACMCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACKCAC", b"GTGCACCCACTGT", 0.0, true, true, true, true, true, false, 1, None);

    //   Reference G
    test_align_both(b"CACBCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACACAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACCCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACDCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACHCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACGCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACTCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACUCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACVCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACRCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACYCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACSCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACWCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACKCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACMCAC", b"GTGCACGCACTGT", 0.0, true, true, true, true, true, false, 1, None);

    //   Reference T
    test_align_both(b"CACBCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACACAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACCCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACDCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACGCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACHCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACVCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACTCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACUCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACYCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACRCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACWCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACSCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, None);
    test_align_both(b"CACKCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, Some((0, 7, 3, 10, 7, 0)));
    test_align_both(b"CACMCAC", b"GTGCACTCACTGT", 0.0, true, true, true, true, true, false, 1, None);
}

#[test]
fn error_rates_align() {
    // Error rates
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACAACAACCACCAAC",
                    0.0, false, true, true, false, false, false, 10,
                            Some((0, 20, 8, 28, 20, 0)));
    // One mismatch at position 7
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACCACCAAC",
                    0.0, false, true, true, false, false, false, 10, None);
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACCACCAAC",
                    0.05, false, true, true, false, false, false, 10,
                    Some((0, 20, 8, 28, 19, 1)));
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACCACCAAC",
                    0.1, false, true, true, false, false, false, 10,
                    Some((0, 20, 8, 28, 19, 1)));
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACCACCAA",
                    0.05, false, true, true, false, false, false, 10, None);
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACCACCAA",
                    0.1, false, true, true, false, false, false, 10,
                    Some((0, 19, 8, 27, 18, 1)));
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACA",
                    0.1, false, true, true, false, false, false, 10,
                    Some((0, 11, 8, 19, 10, 1)));
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGAC",
                    0.1, false, true, true, false, false, false, 10,
                    Some((0, 10, 8, 18, 9, 1)));
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGA",
                    0.1, false, true, true, false, false, false, 10, None);

    // Two mismatches (7 and 13) out of 20 at 0.1 error rate
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACTACCAAC",
                    0.1, false, true, true, false, false, false, 10,
                    Some((0, 20, 8, 28, 18, 2)));
    // Two mismatches (7 and 13) out of 19 at 0.1 error rate, no match
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACTACCAA",
                    0.1, false, true, true, false, false, false, 10, None);

    // Wildcard -- asymmetric reference versus query
    //   1 mismatch in 20 at 0.05 error rate
    test_align(        b"CAACCACAACAACCACCAAC",
               b"GTGTTGTGCAACCACGACAACCACCAAC",
               0.05, false, true, true, false, true, false, 10,
               Some((0, 20, 8, 28, 19, 1)));
    //   1 mismatch in 20 but reference N so 1/19
    test_align(        b"CAACCACAACAACCACNAAC",
               b"GTGTTGTGCAACCACGACAACCACCAAC",
               0.05, false, true, true, false, true, false, 10, None);
    //   1 mismatch in 20 and reference Y (not N) so 1/20
    test_align(        b"CAACCACAACAACCACYAAC",
               b"GTGTTGTGCAACCACGACAACCACCAAC",
               0.05, false, true, true, false, true, false, 10,
               Some((0, 20, 8, 28, 19, 1)));
    //   1 mismatch in 20 and N in query so 1/20
    test_align(        b"CAACCACAACAACCACCAAC",
               b"GTGTTGTGCAACCACGACAACNACCAAC",
               0.05, false, true, true, false, false, true, 10,
               Some((0, 20, 8, 28, 19, 1)));
    //   1 mismatch in 21 and reference N so 1/20, matches
    test_align(        b"CAACCACAACAACCACNAACC",
               b"GTGTTGTGCAACCACGACAACCACCAACC",
               0.05, false, true, true, false, true, false, 10,
               Some((0, 21, 8, 29, 20, 1)));   
}

#[test]
fn edges_align() {
    // Testing edges
    test_align_both(b"AGAGGAG", b"AGGAGTCT",
                    0.0, true, true, true, true, false, false, 1,
                    Some((2, 7, 0, 5, 5, 0)));
    test_align_both(b"AGAGGAG", b"GAGGAGTCT",
                    0.0, true, true, true, true, false, false, 1,
                    Some((1, 7, 0, 6, 6, 0)));
    test_align_both(b"AGAGGAG", b"AGAGGAGTCT",
                    0.0, true, true, true, true, false, false, 1,
                    Some((0, 7, 0, 7, 7, 0)));
    test_align_both(b"AGAGGAG", b"TAGAGGAGTCT",
                    0.0, true, true, true, true, false, false, 1,
                    Some((0, 7, 1, 8, 7, 0)));
    test_align_both(b"AGAGGAG", b"CTAGAGGAGTCT",
                    0.0, true, true, true, true, false, false, 1,
                    Some((0, 7, 2, 9, 7, 0)));
    test_align_both(b"AGAGGAG", b"CTAGAGGAG",
                    0.0, true, true, true, true, false, false, 1,
                    Some((0, 7, 2, 9, 7, 0)));
    test_align_both(b"AGAGGAG", b"TCTAGAGGA",
                    0.0, true, true, true, true, false, false, 1,
                    Some((0, 6, 3, 9, 6, 0)));
    test_align_both(b"AGAGGAG", b"TCTAGAGG",
                    0.0, true, true, true, true, false, false, 1,
                    Some((0, 5, 3, 8, 5, 0)));
}

const INDEL_COST: isize = 1;

fn test_align_both(reference: &[u8],
                   query: &[u8],
                   max_err: f64,
                   start_in_ref: bool,
                   start_in_query: bool,
                   stop_in_ref: bool,
                   stop_in_query: bool,
                   wildcard_ref: bool,
                   wildcard_query: bool,
                   min_overlap: usize,
                   expected: Option<(isize, isize, isize, isize, isize, isize)>)
{
    let exp_flipped = if let Some((refstart, refstop, querystart, querystop, matches, errors)) = expected {
        Some((querystart, querystop, refstart, refstop, matches, errors))
    } else {
        None
    };

    test_align(reference, query, max_err, start_in_ref, start_in_query, stop_in_ref, stop_in_query, wildcard_ref, wildcard_query, min_overlap, expected);
    test_align(query, reference, max_err, start_in_query, start_in_ref, stop_in_query, stop_in_ref, wildcard_query, wildcard_ref, min_overlap, exp_flipped);
}

fn test_align(reference: &[u8],
              query: &[u8],
              max_err: f64,
              start_in_ref: bool,
              start_in_query: bool,
              stop_in_ref: bool,
              stop_in_query: bool,
              wildcard_ref: bool,
              wildcard_query: bool,
              min_overlap: usize,
              expected: Option<(isize, isize, isize, isize, isize, isize)>)
{
    print!("Aligning reference {:?} vs query {:?} error {:?} bounds {:?} wildcards {:?}\n",
           std::str::from_utf8(reference).unwrap(),
           std::str::from_utf8(query).unwrap(),
           max_err,
           (start_in_ref, start_in_query, stop_in_ref, stop_in_query),
           (wildcard_ref, wildcard_query)
           );
    let flags = (if start_in_ref { 1 } else { 0 })
        + (if start_in_query { 2 } else { 0 })
        + (if stop_in_ref { 4 } else { 0 })
        + (if stop_in_query { 8 } else { 0 });
    let mut aligner = Aligner::init(reference, max_err, flags, wildcard_ref, wildcard_query, INDEL_COST, min_overlap as isize);
    let actual = aligner.locate(query);

    if actual != expected {
        aligner.enable_debug();
        aligner.locate(query);
        print!("Expected {:?}\n", expected);
        print!("Actual   {:?}\n", actual);
        print!("{}\n", aligner.dpmatrix().as_ref().unwrap());
        print!("aligner = {:?}\n", aligner);
    }

    let ref_ends = match (start_in_ref, stop_in_ref) {
        (true, true)   => align::AlignEnds::Local,
        (true, false)  => align::AlignEnds::LocalStart,
        (false, true)  => align::AlignEnds::LocalStop,
        (false, false) => align::AlignEnds::Global,
    };

    let query_ends = match (start_in_query, stop_in_query) {
        (true, true)   => align::AlignEnds::Local,
        (true, false)  => align::AlignEnds::LocalStart,
        (false, true)  => align::AlignEnds::LocalStop,
        (false, false) => align::AlignEnds::Global,
    };

    let matching = match (wildcard_ref, wildcard_query) {
        (false, false) => align::AlignMatching::NoWildcard,
        (true, false)  => align::AlignMatching::RefWildcard,
        (false, true)  => align::AlignMatching::QueryWildcard,
        (true, true)   => panic!("wildcards for ref and query!"),
    };

    // locate(query) -> (refstart, refstop, querystart, querystop, matches, errors)
    let new_aligner_conf
        = align::AlignerConf {
            max_error_rate: max_err,
            reference_ends: ref_ends,
            query_ends: query_ends,
            matching: matching,
            indel_cost: INDEL_COST as usize,
            min_overlap: min_overlap,
        };
    let mut new_aligner = align::Aligner::new(&new_aligner_conf, reference).unwrap();
    let new_actual = new_aligner.locate(query);

    let new_actual_equals_expected
        = if let Some(new_actual_location) = &new_actual {
            if let Some((refstart, refstop, querystart, querystop, matches, errors)) = expected {
                refstart == new_actual_location.refstart() as isize
                    && refstop == new_actual_location.refstop() as isize
                    && querystart == new_actual_location.querystart() as isize
                    && querystop == new_actual_location.querystop() as isize
                    && matches == new_actual_location.matches() as isize
                    && errors == new_actual_location.errors() as isize
            } else {
                false
            }
        } else {
            expected.is_none()
        };
    
    if !new_actual_equals_expected {
        new_aligner.enable_debug();
        new_aligner.locate(query);
        print!("Expected {:?}\n", expected);
        print!("Actual   {:?}\n", actual);
        print!("New      {:?}\n", new_actual);
        print!("{}\n", new_aligner.dpmatrix().as_ref().unwrap());
        print!("new_aligner = {:?}\n", new_aligner);
        std::process::exit(1);
    }
}

#[test]
fn sample_align()
{
    let mut aligner = Aligner::init(b"CACCA", 0.0, 0xf, false, false, 1, 1);
    aligner.enable_debug();
    let loc = aligner.locate(b"GACCACCATTA").unwrap();
    print!("{:?}\n", loc);
    print!("{}\n", aligner.dpmatrix().as_ref().unwrap());
    // print!("{:?}", aligner);
    // assert!(loc.refstart() == 0);
    // assert!(loc.refstop() == 5);
    // assert!(loc.querystart() == 3);
    // assert!(loc.querystop() == 8);
    // assert!(loc.matches() == 5);
    // assert!(loc.errors() == 0);
}


                                                                                // # structure for a DP matrix entry
#[derive(Clone,Copy,Debug,Default)]
struct Entry {                                                                  // ctypedef struct _Entry:
  cost: isize,                                                                  //     int cost
  matches: isize,                                                               //     int matches  # no. of matches in this alignment
  origin: isize,                                                                //     int origin   # where the alignment originated: negative for positions within seq1, positive for pos. within seq2
}

#[derive(Clone,Copy,Debug, Default)]
struct Match {                                                                  // ctypedef struct _Match:
  origin: isize,                                                                //     int origin
  cost: isize,                                                                  //     int cost
  matches: isize,                                                               //     int matches
  ref_stop: isize,                                                              //     int ref_stop
  query_stop: isize,                                                            //     int query_stop
}

#[derive(Clone,Debug)]
pub struct DPMatrix<'a> {                                                       // class DPMatrix:
                                                                                //     """
                                                                                //     Representation of the dynamic-programming matrix.

                                                                                //     This is used only when debugging is enabled in the Aligner class since the
                                                                                //     matrix is normally not stored in full.

                                                                                //     Entries in the matrix may be None, in which case that value was not
                                                                                //     computed.
                                                                                //     """

    rows: Vec<Vec<Option<isize>>>,
    reference: &'a [u8],
    query: &'a [u8],
}

impl <'a> DPMatrix<'a> {
    pub fn init(reference: &'a [u8], query: &'a [u8])                           //     def __init__(self, reference, query):
                -> Self {                         
        let m = reference.len();                                                //         m = len(reference)
        let n = query.len();                                                    //         n = len(query)
        let rows = vec![vec![None; n+1]; m+1];                                  //         self._rows = [ [None] * (n+1) for _ in range(m + 1) ]
        Self { rows: rows,
               reference: reference,                                            //         self.reference = reference
               query: query,                                                    //         self.query = query
        }
    }
    
    pub fn set_entry(&mut self, i: isize, j: isize, cost: isize)                //     def set_entry(self, int i, int j, cost):
                                                                                //         """
                                                                                //         Set an entry in the dynamic programming matrix.
                                                                                //         """
    { 
        self.rows[i as usize][j as usize] = Some(cost);                         //         self._rows[i][j] = cost
    }
}

impl <'a> std::fmt::Display for DPMatrix<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>)                              //     def __str__(self):
           -> std::fmt::Result {         
                                                                                //         """
                                                                                //         Return a representation of the matrix as a string.
                                                                                //         """
        let headers = self.query.iter()                                         //         rows = ['     ' + ' '.join(c.rjust(2) for c in self.query)]
            .map(|c| format!("{:>2}", *c as char))
            .collect::<Vec<String>>();
        let header = "     ".to_string()
            + &headers.join(" ");
        let mut rows = vec![header];
        for (q, row) in std::iter::once(b' ')                                   //         for c, row in zip(' ' + self.reference, self._rows):
            .chain(self.reference.iter().copied())
            .zip(self.rows.iter()) {
                let cells = row.iter()                                          //             r = c + ' ' + ' '.join('  ' if v is None else '{:2d}'.format(v) for v in row)
                    .map(|v| match v {
                        None => "  ".to_string(),
                        Some(c) => format!("{:>2}", c),
                    })
                    .collect::<Vec<String>>();
                let r = format!("{} {}",
                                q as char,
                                cells.join(" "));
                
                rows.push(r);                                                   //             rows.append(r)
            }
        write!(f, "{}", rows.join("\n"))                                        //         return '\n'.join(rows)
    }
}

#[derive(Clone, Debug)]
pub struct Aligner<'a> {                                                        // cdef class Aligner:
                                                                                //     """
                                                                                //     Find a full or partial occurrence of a query string in a reference string
                                                                                //     allowing errors (mismatches, insertions, deletions).

                                                                                //     By default, unit costs are used, meaning that mismatches, insertions and
                                                                                //     deletions are counted as one error (edit distance).

                                                                                //     Semi-global alignments allow skipping a suffix and/or prefix of query or
                                                                                //     reference at no cost. Combining semi-global alignment with edit distance is
                                                                                //     a bit unusual because the trivial “optimal” solution at edit distance 0
                                                                                //     would be to skip all of the reference and all of the query, like this:

                                                                                //         REFERENCE-----
                                                                                //         ---------QUERY

                                                                                //     Conceptually, the algorithm used here instead tests all possible overlaps
                                                                                //     between the two sequences and chooses the overlap which maximizes the
                                                                                //     number of matches in the overlapping part, while the error rate must not
                                                                                //     go above a threshold.

                                                                                //     TODO working here

                                                                                //     To allow skipping of a prefix of string1 at no cost, set the
                                                                                //     START_IN_REFERENCE flag.
                                                                                //     To allow skipping of a prefix of string2 at no cost, set the
                                                                                //     START_IN_QUERY flag.
                                                                                //     If both are set, a prefix of string1 or of string1 is skipped,
                                                                                //     never both.
                                                                                //     Similarly, set STOP_IN_REFERENCE and STOP_IN_QUERY to
                                                                                //     allow skipping of suffixes of string1 or string2. Again, when both
                                                                                //     flags are set, never suffixes in both strings are skipped.
                                                                                //     If all flags are set, this results in standard semiglobal alignment.

                                                                                //     The skipped parts are described with two intervals (start1, stop1),
                                                                                //     (start2, stop2).

                                                                                //     For example, an optimal semiglobal alignment of SISSI and MISSISSIPPI looks like this:

                                                                                //     ---SISSI---
                                                                                //     MISSISSIPPI

                                                                                //     start1, stop1 = 0, 5
                                                                                //     start2, stop2 = 3, 8
                                                                                //     (with zero errors)

                                                                                //     The aligned parts are string1[start1:stop1] and string2[start2:stop2].

                                                                                //     The error rate is: errors / length where length is (stop1 - start1).

                                                                                //     An optimal alignment fulfills all of these criteria:

                                                                                //     - its error_rate is at most max_error_rate
                                                                                //     - Among those alignments with error_rate <= max_error_rate, the alignment contains
                                                                                //       a maximal number of matches (there is no alignment with more matches).
                                                                                //     - If there are multiple alignments with the same no. of matches, then one that
                                                                                //       has minimal no. of errors is chosen.
                                                                                //     - If there are still multiple candidates, choose the alignment that starts at the
                                                                                //       leftmost position within the read.

                                                                                //     The alignment itself is not returned, only the tuple
                                                                                //     (start1, stop1, start2, stop2, matches, errors), where the first four fields have the
                                                                                //     meaning as described, matches is the number of matches and errors is the number of
                                                                                //     errors in the alignment.

                                                                                //     It is always the case that at least one of start1 and start2 is zero.

                                                                                //     IUPAC wildcard characters can be allowed in the reference and the query
                                                                                //     by setting the appropriate flags.

                                                                                //     If neither flag is set, the full ASCII alphabet is used for comparison.
                                                                                //     If any of the flags is set, all non-IUPAC characters in the sequences
                                                                                //     compare as 'not equal'.
                                                                                //     """
                                                                                //     cdef:
    m: isize,                                                                   //         int m
    column: Vec<Entry>,                                                         //         _Entry* column  # one column of the DP matrix
    max_error_rate: f64,                                                        //         double max_error_rate
    start_in_reference: bool,                                                   //         bint start_in_reference
    start_in_query: bool,                                                       //         bint start_in_query
    stop_in_reference: bool,                                                    //         bint stop_in_reference
    stop_in_query: bool,                                                        //         bint stop_in_query
    insertion_cost: isize,                                                      //         int _insertion_cost
    deletion_cost: isize,                                                       //         int _deletion_cost
    min_overlap: isize,                                                         //         int _min_overlap
    wildcard_ref: bool,                                                         //         bint wildcard_ref
    wildcard_query: bool,                                                       //         bint wildcard_query
    debug: bool,                                                                //         bint debug
    dpmatrix: Option<DPMatrix<'a>>,                                             //         object _dpmatrix
    reference: &'a [u8],                                                        //         str reference  # reference as set by the user (as str)
    breference: Vec<u8>,                                                        //         bytes _reference  # internal, bytes version of reference (possibly translated to a non-ASCII representation)
    effective_length: isize,                                                    //         readonly int effective_length
    n_counts: Vec<isize>,                                                       //         int* n_counts  # n_counts[i] == number of N characters in reference[:i]
}

impl <'a> Aligner<'a> {
    pub fn init(                                                                //     def __cinit__(
                                                                                //         self,
        reference: &'a [u8],                                                    //         str reference,
        max_error_rate: f64,                                                    //         double max_error_rate,
        flags: isize,                                                           //         int flags=15,
        wildcard_ref: bool,                                                     //         bint wildcard_ref=False,
        wildcard_query: bool,                                                   //         bint wildcard_query=False,
        indel_cost: isize,                                                      //         int indel_cost=1,
        min_overlap: isize,                                                     //         int min_overlap=1,
    ) -> Self                                                                   //     ):
    {
        let mut aligner = Aligner {
            m: 0,
            column: Vec::new(),
            max_error_rate: max_error_rate,                                     //         self.max_error_rate = max_error_rate
            start_in_reference: flags & 1 > 0,                                  //         self.start_in_reference = flags & 1
            start_in_query: flags & 2 > 0,                                      //         self.start_in_query = flags & 2
            stop_in_reference: flags & 4 > 0,                                   //         self.stop_in_reference = flags & 4
            stop_in_query: flags & 8 > 0,                                       //         self.stop_in_query = flags & 8
            wildcard_ref: wildcard_ref,                                         //         self.wildcard_ref = wildcard_ref
            wildcard_query: wildcard_query,                                     //         self.wildcard_query = wildcard_query
                                                                                // ZZZ //         self._set_reference(reference)
                                                                                //         if min_overlap < 1:
                                                                                //             raise ValueError('min_overlap must be at least 1')
            min_overlap: min_overlap,                                           //         self._min_overlap = min_overlap
            debug: false,                                                       //         self.debug = False
            dpmatrix: None,                                                     //         self._dpmatrix = None
            reference: reference,
            breference: Vec::with_capacity(reference.len()),
            effective_length: 0,
            n_counts: Vec::new(),
                                                                                //         if indel_cost < 1:
                                                                                //             raise ValueError('indel_cost must be at least 1')
            insertion_cost: indel_cost,                                         //         self._insertion_cost = indel_cost
            deletion_cost: indel_cost, };                                       //         self._deletion_cost = indel_cost
        aligner.set_reference(reference);
        aligner
    }
    
    fn set_reference(&mut self, reference: &'a [u8])                            //     def _set_reference(self, str reference):
    {

        let mem: Vec<Entry> = vec![Default::default(); reference.len() + 1];    //         mem = <_Entry*> PyMem_Realloc(self.column, (len(reference) + 1) * sizeof(_Entry))
                                                                                //         if not mem:
                                                                                //             raise MemoryError()
        let mem_nc: Vec<isize> = vec![Default::default(); reference.len() + 1]; //         mem_nc = <int*> PyMem_Realloc(self.n_counts, (len(reference) + 1) * sizeof(int))
                                                                                //         if not mem_nc:
                                                                                //             raise MemoryError()
        self.column = mem;                                                      //         self.column = mem
        self.n_counts = mem_nc;                                                 //         self.n_counts = mem_nc
        self.breference = reference.to_vec();                                   //         self._reference = reference.encode('ascii')
        self.m = reference.len() as isize;                                      //         self.m = len(reference)
        self.effective_length = self.m;                                         //         self.effective_length = self.m
        let mut n_count = 0;                                                    //         n_count = 0
        for i in 0..self.m {                                                    //         for i in range(self.m):
            self.n_counts[i as usize] = n_count;                                //             self.n_counts[i] = n_count
            if reference[i as usize] == b'n' || reference[i as usize] == b'N' { //             if reference[i] == 'n' or reference[i] == 'N':
                n_count += 1;                                                   //                 n_count += 1
            }
        }
        self.n_counts[self.m as usize] = n_count;                               //         self.n_counts[self.m] = n_count
        assert_eq!(self.n_counts[self.m as usize],                              //         assert self.n_counts[self.m] == reference.count('N') + reference.count('n')
                   reference.iter().copied()
                   .filter(|&c| c == b'N' || c == b'n').count() as isize); 
        if self.wildcard_ref {                                                  //         if self.wildcard_ref:
            self.effective_length = self.m - self.n_counts[self.m as usize];    //             self.effective_length = self.m - self.n_counts[self.m]
            if self.effective_length == 0 {                                     //             if self.effective_length == 0:
                panic!("Cannot have only N wildcards in the sequence");         //                 raise ValueError("Cannot have only N wildcards in the sequence")
            }
            encode_iupac_vec(&self.reference, &mut self.breference);            //             self._reference = self._reference.translate(IUPAC_TABLE)
        } else if self.wildcard_query {                                         //         elif self.wildcard_query:
            encode_acgt_vec(&self.reference, &mut self.breference);             //             self._reference = self._reference.translate(ACGT_TABLE)
        }
        self.reference = reference;                                         //         self.reference = reference
    }   
                                                                                //     property dpmatrix:
                                                                                //         """
                                                                                //         The dynamic programming matrix as a DPMatrix object. This attribute is
                                                                                //         usually None, unless debugging has been enabled with enable_debug().
                                                                                //         """
    pub fn dpmatrix(&self) -> &Option<DPMatrix> {                               //         def __get__(self):
        &self.dpmatrix                                                          //             return self._dpmatrix
    }
    
    pub fn enable_debug(&mut self) {                                            //     def enable_debug(self):
                                                                                //         """
                                                                                //         Store the dynamic programming matrix while running the locate() method
                                                                                //         and make it available in the .dpmatrix attribute.
                                                                                //         """
        self.debug = true;                                                      //         self.debug = True
    }



    pub fn locate(&mut self, query: &'a [u8])                                   //     def locate(self, str query):
                  -> Option<(isize, isize, isize, isize, isize, isize)> {
                                                                                //         """
                                                                                //         locate(query) -> (refstart, refstop, querystart, querystop, matches, errors)

                                                                                //         Find the query within the reference associated with this aligner. The
                                                                                //         intervals (querystart, querystop) and (refstart, refstop) give the
                                                                                //         location of the match.

                                                                                //         That is, the substrings query[querystart:querystop] and
                                                                                //         self.reference[refstart:refstop] were found to align best to each other,
                                                                                //         with the given number of matches and the given number of errors.

                                                                                //         The alignment itself is not returned.
                                                                                //         """
                                                                                //         cdef:
        let s1 = &self.breference;                                                //             char* s1 = self._reference
        let mut query_bytes = query.to_vec();                                   //             bytes query_bytes = query.encode('ascii')
        let s2: &[u8];                                                          //             char* s2
        let m = self.m;                                                         //             int m = self.m
        let n = query.len() as isize;                                           //             int n = len(query)
        let column = &mut self.column;                                          //             _Entry* column = self.column  # Current column of the DP matrix
        let max_error_rate = self.max_error_rate;                               //             double max_error_rate = self.max_error_rate
        let stop_in_query = self.stop_in_query;                                 //             bint stop_in_query = self.stop_in_query
        let mut compare_ascii = false;                                          //             bint compare_ascii = False
        
        if self.wildcard_query {                                                //         if self.wildcard_query:
            encode_iupac_vec(query, &mut query_bytes);                          //             query_bytes = query_bytes.translate(IUPAC_TABLE)
        } else if self.wildcard_ref {                                           //         elif self.wildcard_ref:
            encode_acgt_vec(query, &mut query_bytes);                           //             query_bytes = query_bytes.translate(ACGT_TABLE)
        } else {                                                                //         else:
                                                                                //             # TODO Adding the .upper() increases overall runtime slightly even
                                                                                //             # when I remove the .upper() from Adapter.match_to().
                                                                                //             query_bytes = query_bytes.upper()
            compare_ascii = true;                                               //             compare_ascii = True
        }
        s2 = &query_bytes;
                                                                                //         s2 = query_bytes
                                                                                //         """
                                                                                //         DP Matrix:
                                                                                //                    query (j)
                                                                                //                  ----------> n
                                                                                //                 |
                                                                                //         ref (i) |
                                                                                //                 |
                                                                                //                 V
                                                                                //                m
                                                                                //         """
        
                                                                                //         # maximum no. of errors
        let k = (max_error_rate * m as f64).floor() as isize;                   //         cdef int k = <int> (max_error_rate * m)

                                                                                //         # Determine largest and smallest column we need to compute
        let mut max_n = n;                                                      //         cdef int max_n = n
        let mut min_n = 0;                                                      //         cdef int min_n = 0
        if !self.start_in_query {                                               //         if not self.start_in_query:
                                                                                //             # costs can only get worse after column m
            max_n = isize::min(n, m + k);                                       //             max_n = min(n, m + k)
        }
        if !self.stop_in_query {                                                //         if not self.stop_in_query:
            min_n = isize::max(0, n - m - k);                                   //             min_n = max(0, n - m - k)
        }
        
                                                                                //         # Fill column min_n.
                                                                                //         #
                                                                                //         # Four cases:
                                                                                //         # not startin1, not startin2: c(i,j) = max(i,j); origin(i, j) = 0
                                                                                //         #     startin1, not startin2: c(i,j) = j       ; origin(i, j) = min(0, j - i)
                                                                                //         # not startin1,     startin2: c(i,j) = i       ; origin(i, j) =
                                                                                //         #     startin1,     startin2: c(i,j) = min(i,j)

                                                                                //         # TODO (later)
                                                                                //         # fill out columns only until 'last'
        if !self.start_in_reference && !self.start_in_query {                   //         if not self.start_in_reference and not self.start_in_query:
            for i in 0..(m+1) {                                                 //             for i in range(m + 1):
                column[i as usize].matches = 0;                                 //                 column[i].matches = 0
                column[i as usize].cost =                                       //                 column[i].cost = max(i, min_n) * self._insertion_cost
                    isize::max(i, min_n) * self.insertion_cost; 
                column[i as usize].origin = 0;                                  //                 column[i].origin = 0
            }
        } else if self.start_in_reference && !self.start_in_query {             //         elif self.start_in_reference and not self.start_in_query:
            for i in 0..(m+1) {                                                 //             for i in range(m + 1):
                column[i as usize].matches = 0;                                 //                 column[i].matches = 0
                column[i as usize].cost = min_n * self.insertion_cost;          //                 column[i].cost = min_n * self._insertion_cost
                column[i as usize].origin = isize::min(0, min_n - i);           //                 column[i].origin = min(0, min_n - i)
            }
        } else if !self.start_in_reference && self.start_in_query {             //         elif not self.start_in_reference and self.start_in_query:
            for i in 0..(m+1) {                                                 //             for i in range(m + 1):
                column[i as usize].matches = 0;                                 //                 column[i].matches = 0
                column[i as usize].cost = i * self.insertion_cost;              //                 column[i].cost = i * self._insertion_cost
                column[i as usize].origin = isize::max(0, min_n - i);           //                 column[i].origin = max(0, min_n - i)
            }
        } else {                                                                //         else:
            for i in 0..(m+1) {                                                 //             for i in range(m + 1):
                column[i as usize].matches = 0;                                 //                 column[i].matches = 0
                column[i as usize].cost                                         //                 column[i].cost = min(i, min_n) * self._insertion_cost
                    = isize::min(i, min_n) * self.insertion_cost; 
                column[i as usize].origin = min_n - i;                          //                 column[i].origin = min_n - i
            }
        }

        if self.debug {                                                         //         if self.debug:
            let mut dpmatrix = DPMatrix::init(self.reference, query);           //             self._dpmatrix = DPMatrix(self.reference, query)
            for i in 0..(m+1) {                                                 //             for i in range(m + 1):
                dpmatrix.set_entry(i, min_n, column[i as usize].cost);          //                 self._dpmatrix.set_entry(i, min_n, column[i].cost)
            }
            self.dpmatrix = Some(dpmatrix);
        }
        
        let mut best = Match::default();                                        //         cdef _Match best
        best.ref_stop = m;                                                      //         best.ref_stop = m
        best.query_stop = n;                                                    //         best.query_stop = n
        best.cost = m + n;                                                      //         best.cost = m + n
        best.origin = 0;                                                        //         best.origin = 0
        best.matches = 0;                                                       //         best.matches = 0

                                                                                //         # Ukkonen's trick: index of the last cell that is at most k
        let mut last = isize::min(m, k + 1);                                    //         cdef int last = min(m, k + 1)
        if self.start_in_reference {                                            //         if self.start_in_reference:
            last = m;                                                           //             last = m
        }
        
                                                                                //         cdef:
        let mut cost_diag;                                                      //             int cost_diag
        let mut cost_deletion;                                                  //             int cost_deletion
        let mut cost_insertion;                                                 //             int cost_insertion
        let mut origin;                                                         //             int origin, cost, matches
        let mut cost;
        let mut matches;
        let mut length;                                                         //             int length
        let mut ref_start: isize;                                               //             int ref_start
        let mut cur_effective_length;                                           //             int cur_effective_length
        let mut characters_equal;                                               //             bint characters_equal
                                                                                //             # We keep only a single column of the DP matrix in memory.
                                                                                //             # To access the diagonal cell to the upper left,
                                                                                //             # we store it here before overwriting it.
        let mut diag_entry;                                                     //             _Entry diag_entry

        {                                                                       //         with nogil:
                                                                                //             # iterate over columns
            for j in (min_n+1)..(max_n+1) {                                     //             for j in range(min_n + 1, max_n + 1):
                                                                                //                 # remember first entry before overwriting
                diag_entry = column[0];                                         //                 diag_entry = column[0]
            
        
                                                                                //                 # fill in first entry in this column
                if self.start_in_query {                                        //                 if self.start_in_query:
                    column[0].origin = j;                                       //                     column[0].origin = j
                } else {                                                        //                 else:
                    column[0].cost = j * self.insertion_cost;                   //                     column[0].cost = j * self._insertion_cost
                }

                for i in 1..(last+1) {                                          //                 for i in range(1, last + 1):
                    if compare_ascii {                                          //                     if compare_ascii:
                        characters_equal =                                      //                         characters_equal = (s1[i-1] == s2[j-1])
                            s1[(i-1) as usize] == s2[(j-1) as usize]; 
                    } else {                                                    //                     else:
                        characters_equal =                                      //                         characters_equal = (s1[i-1] & s2[j-1]) != 0
                            (s1[(i-1) as usize] & s2[(j-1) as usize]) != 0; 
                    }
                    
                    if characters_equal {                                       //                     if characters_equal:
                                                                                //                         # If the characters match, skip computing costs for
                                                                                //                         # insertion and deletion as they are at least as high.
                        cost = diag_entry.cost;                                 //                         cost = diag_entry.cost
                        origin = diag_entry.origin;                             //                         origin = diag_entry.origin
                        matches = diag_entry.matches + 1;                       //                         matches = diag_entry.matches + 1
                    } else {                                                    //                     else:
                                                                                //                         # Characters do not match.
                        cost_diag = diag_entry.cost + 1;                        //                         cost_diag = diag_entry.cost + 1
                        cost_deletion =                                         //                         cost_deletion = column[i].cost + self._deletion_cost
                            column[i as usize].cost + self.deletion_cost;
                        cost_insertion =                                        //                         cost_insertion = column[i-1].cost + self._insertion_cost
                            column[(i-1) as usize].cost + self.insertion_cost; 

                        if cost_diag <= cost_deletion                           //                         if cost_diag <= cost_deletion and cost_diag <= cost_insertion:
                            && cost_diag <= cost_insertion {
                                                                                //                             # MISMATCH
                                cost = cost_diag;                               //                             cost = cost_diag
                                origin = diag_entry.origin;                     //                             origin = diag_entry.origin
                                matches = diag_entry.matches;                   //                             matches = diag_entry.matches
                            } else if cost_insertion < cost_deletion {          //                         elif cost_insertion <= cost_deletion:
                                                                                //                             # INSERTION
                                cost = cost_insertion;                          //                             cost = cost_insertion
                                origin = column[(i-1) as usize].origin;         //                             origin = column[i-1].origin
                                matches = column[(i-1) as usize].matches;       //                             matches = column[i-1].matches
                            } else {                                            //                         else:
                                                                                //                             # DELETION
                                cost = cost_deletion;                           //                             cost = cost_deletion
                                origin = column[i as usize].origin;             //                             origin = column[i].origin
                                matches = column[i as usize].matches;           //                             matches = column[i].matches
                            }
                    }
                                                                                //                     # Remember the current cell for next iteration
                    diag_entry = column[i as usize];                            //                     diag_entry = column[i]
                    
                    column[i as usize].cost = cost;                             //                     column[i].cost = cost
                    column[i as usize].origin = origin;                         //                     column[i].origin = origin
                    column[i as usize].matches = matches;                       //                     column[i].matches = matches
                }
                if self.debug {                                                 //                 if self.debug:
                    {                                                           //                     with gil:
                        for i in 0..(last+1) {                                  //                         for i in range(last + 1):
                            self.dpmatrix                                       //                             self._dpmatrix.set_entry(i, j, column[i].cost)
                                .as_mut()
                                .map(|dpm| dpm.set_entry(i, j, column[i as usize].cost));
                        }
                    }
                }
                while last >= 0 && column[last as usize].cost > k {             //                 while last >= 0 and column[last].cost > k:
                    last -= 1;                                                  //                     last -= 1
                }
            
                                                                                //                 # last can be -1 here, but will be incremented next.
                                                                                //                 # TODO if last is -1, can we stop searching?
                if last < m {                                                   //                 if last < m:
                    last += 1;                                                  //                     last += 1
                } else if stop_in_query {                                       //                 elif stop_in_query:
                                                                                //                     # Found a match. If requested, find best match in last row.
                                                                                //                     # length of the aligned part of the reference
                    length = m + isize::min(column[m as usize].origin, 0);      //                     length = m + min(column[m].origin, 0)
                    cur_effective_length = length;                              //                     cur_effective_length = length
                    if self.wildcard_ref {                                      //                     if self.wildcard_ref:
                        if length < m {                                         //                         if length < m:
                                                                                //                             # Recompute effective length so that it only takes into
                                                                                //                             # account the matching suffix of the reference
                            cur_effective_length =                              //                             cur_effective_length = length - self.n_counts[length]
                                length - self.n_counts[length as usize];
                        } else {                                                //                         else:
                            cur_effective_length = self.effective_length;       //                             cur_effective_length = self.effective_length
                        }
                    }
                    cost = column[m as usize].cost;                             //                     cost = column[m].cost
                    matches = column[m as usize].matches;                       //                     matches = column[m].matches
                    if length >= self.min_overlap                               //                     if length >= self._min_overlap and cost <= cur_effective_length * max_error_rate and (matches > best.matches or (matches == best.matches and cost < best.cost)):
                        && (cost as f64) <= (cur_effective_length as f64) * max_error_rate
                        && (matches > best.matches
                            || (matches == best.matches && cost < best.cost)) {
                                                                                //                         # update
                            best.matches = matches;                             //                         best.matches = matches
                            best.cost = cost;                                   //                         best.cost = cost
                            best.origin = column[m as usize].origin;            //                         best.origin = column[m].origin
                            best.ref_stop = m;                                  //                         best.ref_stop = m
                            best.query_stop = j;                                //                         best.query_stop = j
                            if cost == 0 && matches == m {                      //                         if cost == 0 and matches == m:
                                                                                //                             # exact match, stop early
                                break;                                          //                             break
                            }
                                                                                //                 # column finished
                        }
                }
            }
        }
        
    
        if max_n == n {                                                         //         if max_n == n:
            let first_i = if self.stop_in_reference { 0 } else { m };           //             first_i = 0 if self.stop_in_reference else m
                                                                                //             # search in last column # TODO last?
            for i in first_i..(m+1) {                                           //             for i in range(first_i, m+1):
                length = i + isize::min(column[i as usize].origin, 0);          //                 length = i + min(column[i].origin, 0)
                cost = column[i as usize].cost;                                 //                 cost = column[i].cost
                matches = column[i as usize].matches;                           //                 matches = column[i].matches
                if self.wildcard_ref {                                          //                 if self.wildcard_ref:
                    if length < m {                                             //                     if length < m:
                                                                                //                         # Recompute effective length so that it only takes into
                                                                                //                         # account the matching part of the reference
                        ref_start =                                             //                         ref_start = -min(column[i].origin, 0)
                            - isize::min(column[i as usize].origin,0);
                        assert!(0 <= ref_start);                                //                         assert 0 <= ref_start <= m
                        assert!(ref_start <= m);
                        cur_effective_length =                                  //                         cur_effective_length = length - (self.n_counts[i] - self.n_counts[ref_start])
                            length - (self.n_counts[i as usize] - self.n_counts[ref_start as usize]);
                    } else {                                                    //                     else:
                        cur_effective_length = self.effective_length;           //                         cur_effective_length = self.effective_length
                    }
                } else {                                                        //                 else:
                    cur_effective_length = length;                              //                     cur_effective_length = length
                }
            
                assert!(0 <= cur_effective_length);                             //                 assert 0 <= cur_effective_length and cur_effective_length <= length
                assert!(cur_effective_length <= length);
                assert!(cur_effective_length <= self.effective_length);         //                 assert cur_effective_length <= self.effective_length

                if length >= self.min_overlap                                   //                 if length >= self._min_overlap and cost <= cur_effective_length * max_error_rate and (matches > best.matches or (matches == best.matches and cost < best.cost)):
                    && (cost as f64) <= cur_effective_length as f64 * max_error_rate
                    && (matches > best.matches
                        || (matches == best.matches && cost < best.cost)) {
                                                                                //                     # update best
                        best.matches = matches;                                 //                     best.matches = matches
                        best.cost = cost;                                       //                     best.cost = cost
                        best.origin = column[i as usize].origin;                //                     best.origin = column[i].origin
                        best.ref_stop = i;                                      //                     best.ref_stop = i
                        best.query_stop = n;                                    //                     best.query_stop = n
                    }
            }
        }
        
        if best.cost == m + n {                                                 //         if best.cost == m + n:
                                                                                //             # best.cost was initialized with this value.
                                                                                //             # If it is unchanged, no alignment was found that has
                                                                                //             # an error rate within the allowed range.
            return None;                                                        //             return None
        }
        
        let start1;                                                             //         cdef int start1, start2
        let start2;
        if best.origin >= 0 {                                                   //         if best.origin >= 0:
            start1 = 0;                                                         //             start1 = 0
            start2 = best.origin;                                               //             start2 = best.origin
        } else {                                                                //         else:
            start1 = -best.origin;                                              //             start1 = -best.origin
            start2 = 0;                                                         //             start2 = 0
        }
        assert!(best.ref_stop - start1 > 0);                                    //         assert best.ref_stop - start1 > 0  # Do not return empty alignments.
        return Some((start1, best.ref_stop, start2, best.query_stop,            //         return (start1, best.ref_stop, start2, best.query_stop, best.matches, best.cost)
                     best.matches, best.cost)); 
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

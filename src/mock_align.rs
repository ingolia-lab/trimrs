use trimrs::align;
use trimrs::cutadapt_align::*;

// locate(query) -> (refstart, refstop, querystart, querystop, matches, errors)

fn main() {
    sample_align();
    test_align_both(b"CACCA", b"GACCACCATTA",
                    0.00, true, true, true, true, false, false, 1,
                    Some((0, 5, 3, 8, 5, 0)));

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

    // 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1
    
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

    // Two mismatches at position 7 and 13
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACTACCAAC",
                    0.1, false, true, true, false, false, false, 10,
                    Some((0, 20, 8, 28, 18, 2)));
    test_align_both(        b"CAACCACAACAACCACCAAC",
                    b"GTGTTGTGCAACCACGACAACTACCAA",
                    0.1, false, true, true, false, false, false, 10, None);
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

    let new_aligner_conf
        = align::AlignerConf {
            max_error_rate: max_err,
            reference_ends: ref_ends,
            query_ends: query_ends,
            matching: matching,
            indel_cost: INDEL_COST,
            min_overlap: min_overlap as isize,
        };
    let mut new_aligner = align::Aligner::new(&new_aligner_conf, reference).unwrap();
    let new_actual = new_aligner.locate(query);

    let new_actual_equals_expected
        = if let Some(new_actual_location) = &new_actual {
            if let Some((refstart, refstop, querystart, querystop, matches, errors)) = expected {
                refstart == new_actual_location.refstart()
                    && refstop == new_actual_location.refstop()
                    && querystart == new_actual_location.querystart()
                    && querystop == new_actual_location.querystop()
                    && matches == new_actual_location.matches()
                    && errors == new_actual_location.errors()
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

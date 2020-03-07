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
    } else {
        print!("Good {:?} vs {:?} => {:?}\n",
               std::str::from_utf8(query).unwrap(),
               std::str::from_utf8(reference).unwrap(),
               actual);
    }

    let mut new_aligner = align::Aligner::new(reference, max_err, flags, wildcard_ref, wildcard_query, INDEL_COST, min_overlap as isize).unwrap();
    let new_actual = new_aligner.locate(query);
    if new_actual != expected {
        new_aligner.enable_debug();
        new_aligner.locate(query);
        print!("Expected {:?}\n", expected);
        print!("Actual   {:?}\n", actual);
        print!("New      {:?}\n", new_actual);
        print!("{}\n", new_aligner.dpmatrix().as_ref().unwrap());
    } else {
        print!("No change\n");
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

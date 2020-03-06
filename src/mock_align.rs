use trimrs::cutadapt_align::*;

fn main() {
    let mut aligner = Aligner::init(b"CACCA", 0.0, 0xf, false, false, 1, 1);
    aligner.enable_debug();
    let loc = aligner.locate(b"GACCACCATTA").unwrap();
    print!("{:?}", loc);
    print!("{}", aligner.dpmatrix().as_ref().unwrap());
    print!("{:?}", aligner);
    // assert!(loc.refstart() == 0);
    // assert!(loc.refstop() == 5);
    // assert!(loc.querystart() == 3);
    // assert!(loc.querystop() == 8);
    // assert!(loc.matches() == 5);
    // assert!(loc.errors() == 0);
}

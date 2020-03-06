use std::env;
use std::fs;
use std::io;
use std::io::{Read,Write};
use std::path::{PathBuf};

fn main() -> io::Result<()> {
    let args = env::args().collect::<Vec<String>>();
    if args.len() != 2 {
        eprintln!("Exactly one argument, the original filename");
        return Ok(());
    }

    let orig_filename = &args[1];
    let mut bak_filename = PathBuf::from(orig_filename);
    bak_filename.set_extension("bak");
    fs::rename(&orig_filename, &bak_filename)?;
    
    let orig = fs::read_to_string(&bak_filename)?;

    let mut out = fs::File::create(&orig_filename)?;
    
    for orig_line in orig.lines() {
        if let Some(comment_start) = orig_line.find("//") {
            let (before, after) = orig_line.split_at(comment_start);
            write!(out, "{:80}{}\n", before.trim_end(), after)?;
        } else {
            write!(out, "{}\n", orig_line)?;
        }
    }
    
    Ok(())
}

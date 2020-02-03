use std::sync::Arc;

use crate::output::Fate;

#[derive(Clone, Debug)]
pub struct Trimming<'a> {
    name: &'a [u8],
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
    trim_start: usize,
    trim_len: usize,
    umi_indices: Vec<usize>,
    read_tag: Option<Arc<String>>,
    fate: Fate,
}

impl <'a> Trimming<'a> {
    pub fn new(name: &'a [u8], seq: &'a [u8], qual: Option<&'a [u8]>) -> Self {
        Trimming { name: name,
                   seq: seq,
                   qual: qual,
                   trim_start: 0,
                   trim_len: seq.len(),
                   umi_indices: Vec::new(),
                   read_tag: None,
                   fate: Fate::Output,
        }
    }

    #[inline(always)]
    fn trim_end(&self) -> usize { self.trim_start + self.trim_len }
    
    pub fn seq_curr(&self) -> &'a [u8] { &self.seq[self.trim_start..self.trim_end()] }

    pub fn qual_curr(&self) -> Option<&'a [u8]> { self.qual.map(|q| &q[self.trim_start..self.trim_end()]) }

    pub fn len_curr(&self) -> usize { self.trim_len }
    
    pub fn seq_raw(&self) -> &'a [u8] { self.seq }

    pub fn qual_raw(&self) -> Option<&'a [u8]> { self.qual }

    pub fn name_raw(&self) -> &'a [u8] { self.name }

    pub fn trim_from_start(&mut self, len: usize) -> usize {
        let reallen = len.min(self.trim_len);
        self.trim_start += reallen;
        self.trim_len -= reallen;
        reallen
    }

    pub fn trim_from_end(&mut self, len: usize) -> usize {
        let reallen = len.min(self.trim_len);
        self.trim_len -= reallen;
        reallen
    }
}

/// Encode `nts` without IUPAC wildcards.
///
/// # Arguments
/// * `nts` is the ASCII representation of the sequence
/// * `mask` is modified to contain the encoded mask
pub fn encode_acgt_vec(nts: &[u8], mask: &mut Vec<u8>) {
    mask.resize(nts.len(), 0);
    for i in 0..nts.len() {
        mask[i] = encode_acgt(nts[i]);
    }
}

/// Encodes `nts` with IUPAC wildcards.
///
/// Bitwise and of a wildcard-encoded sequence with an ACGT-encoded
/// sequence will give a non-zero result exactly when the wildcard
/// matches the sequence.
///
/// # Arguments
/// * `nts` is the ASCII representation of the sequence with IUPAC wildcards.
/// * `mask` is modified to contain the encoded mask
pub fn encode_iupac_vec(nts: &[u8], mask: &mut Vec<u8>) {
    mask.resize(nts.len(), 0);
    for i in 0..nts.len() {
        mask[i] = encode_iupac(nts[i]);
    }
}

/// Decodes `mask`, using IUPAC wildcards when needed.
///
/// Any `0` values are decoded to `x`.
///
/// # Arguments
/// * `mask` is an encoded mask from [`encode_acgt_vec`](fn.encode_acgt_vec.html) or [`encode_iupac_vec`](fn.encode_iupac_vec.html)
/// * `nts` is the decoded ASCII string.
pub fn decode_iupac_vec(mask: &[u8], nts: &mut Vec<u8>) {
    nts.resize(mask.len(), 0);
    for i in 0..mask.len() {
        nts[i] = decode_iupac(mask[i]);
    }
}

#[inline]
fn encode_acgt(nt: u8) -> u8 {
    NT_ACGT[nt as usize]
}

#[inline]
fn encode_iupac(nt: u8) -> u8 {
    NT_IUPAC[nt as usize]
}

#[inline]
fn decode_iupac(mask: u8) -> u8 {
    NT_DECODE[(mask & DECODE_MASK) as usize]
}

const ENCODE_SIZE: usize = 0x100;
const DECODE_SIZE: usize = 0x10;
const DECODE_MASK: u8 = 0x0f;

const NT_0: u8 = 0x00;
const NT_A: u8 = 0x01;
const NT_C: u8 = 0x02;
const NT_G: u8 = 0x04;
const NT_T: u8 = 0x08;

static NT_ACGT: [u8; ENCODE_SIZE] =
    [ NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x10
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x20
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x30
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x40
      NT_0, NT_A, NT_0, NT_C, NT_0, NT_0, NT_0, NT_G,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x50
      NT_0, NT_0, NT_0, NT_0, NT_T, NT_T, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x60
      NT_0, NT_A, NT_0, NT_C, NT_0, NT_0, NT_0, NT_G,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x70
      NT_0, NT_0, NT_0, NT_0, NT_T, NT_T, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x80
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
    ];

static NT_IUPAC: [u8; ENCODE_SIZE] =
    [ NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x10
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x20
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x30
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x40
      NT_0,
      NT_A,
      NT_C | NT_G | NT_T,
      NT_C,
      NT_A | NT_G | NT_T,
      NT_0,
      NT_0,
      NT_G,
      NT_A | NT_C | NT_T,
      NT_0,
      NT_0,
      NT_G | NT_T,
      NT_0,
      NT_A | NT_C,
      NT_A | NT_C | NT_G | NT_T,
      NT_0,
// 0x50
      NT_0,
      NT_0,
      NT_A | NT_G,
      NT_C | NT_G,
      NT_T,
      NT_T,
      NT_A | NT_C | NT_G,
      NT_A | NT_T,
      NT_0,
      NT_C | NT_T,
      NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x60
      NT_0,
      NT_A,
      NT_C | NT_G | NT_T,
      NT_C,
      NT_A | NT_G | NT_T,
      NT_0,
      NT_0,
      NT_G,
      NT_A | NT_C | NT_T,
      NT_0,
      NT_0,
      NT_G | NT_T,
      NT_0,
      NT_A | NT_C,
      NT_A | NT_C | NT_G | NT_T,
      NT_0,
// 0x70
      NT_0,
      NT_0,
      NT_A | NT_G,
      NT_C | NT_G,
      NT_T,
      NT_T,
      NT_A | NT_C | NT_G,
      NT_A | NT_T,
      NT_0,
      NT_C | NT_T,
      NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0,
// 0x80
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
      NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0, NT_0,
    ];

static NT_DECODE: [u8; DECODE_SIZE] =
    [ b'x',
      b'A', // ...A
      b'C', // ..C.
      b'M', // ..CA
      b'G', // .G..
      b'R', // .G.A
      b'S', // .GC.
      b'V', // .GCA 
      b'T', // T...
      b'W', // T..A
      b'Y', // T.C.
      b'H', // T.CA
      b'K', // TG..
      b'D', // TG.A
      b'B', // TGC.
      b'N', // TGCA
    ];

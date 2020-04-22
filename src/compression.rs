use vbyte::{VByteEncoded, VByteDecoder};
//use bitpacking::{BitPacker4x, BitPacker};
use deflate::deflate_bytes;

pub enum IntegerEncode {
    Uncoded(Vec<u64>),
    Delta(Vec<u64>),
    VByte(Vec<u8>),
    Pfor(Vec<u32>)
}

enum FloatEncode {
    Uncoded(Vec<f64>)
}

pub enum StringEncode {
    Uncoded(Vec<String>),
    Deflate(Vec<String>)
}

pub fn integer_encode(input: Vec<u64>) -> IntegerEncode {
    let mut vec = Vec::with_capacity(input.len()-1);
    for i in 0..input.len()-1 {
        // vec.push(input[i+1] - input[i]);
        VByteEncoded::new(input[i+1] - input[i]).write_to(&mut vec);
    }
    IntegerEncode::VByte(vec)
}

pub fn integer_decode(input: IntegerEncode) -> Vec<u64> {
    match input {
        IntegerEncode::VByte(vec) => VByteDecoder::new(vec.as_slice()).collect::<Vec<_>>(),
        IntegerEncode::Uncoded(vec) => vec,
        _ => panic!("Not implemented!")
    }
}

pub fn string_encode(input: Vec<String>) -> StringEncode {
    
}
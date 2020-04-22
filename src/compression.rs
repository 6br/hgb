use std::io::{self, Read};
use vbyte::{VByteEncoded, VByteDecoder};
//use bitpacking::{BitPacker4x, BitPacker};
use libflate::deflate::{Encoder, Decoder};


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
    Deflate(Vec<u8>)
}

/// Returns an encoded array of u64.
pub fn integer_encode(input: Vec<u64>) -> IntegerEncode {
    let mut vec = Vec::with_capacity(input.len());
    VByteEncoded::new(input[0]).write_to(&mut vec);

    // Delta-encoding
    for i in 0..input.len()-1 {
        // vec.push(input[i+1] - input[i]);
        VByteEncoded::new(input[i+1] - input[i]).write_to(&mut vec);
    }
    IntegerEncode::VByte(vec)
}

/// Returns a decoded array of u64
pub fn integer_decode(input: IntegerEncode) -> Vec<u64> {
    match input {
        IntegerEncode::VByte(vec) => {
            let decoded = VByteDecoder::new(vec.as_slice()).collect::<Vec<_>>();

            // Delta-decoding
            let mut vec = Vec::with_capacity(decoded.len());
            let mut previous_value = decoded[0];
            vec.push(previous_value);
            
            for i in 1..decoded.len() {
                vec.push(decoded[i] + previous_value);
                previous_value = decoded[i];
            }
            vec
        },
        IntegerEncode::Uncoded(vec) => vec,
        _ => panic!("Not implemented!")
    }
}

pub fn string_encode(input: Vec<String>) -> StringEncode {
    let mut concatenate_string = input.join("\0");
    let mut encoder = Encoder::new(Vec::new());
    io::copy(&mut concatenate_string.as_bytes(), &mut encoder).unwrap();
    let encoded_data = encoder.finish().into_result().unwrap();
    StringEncode::Deflate(encoded_data)
}

pub fn string_decode(input: StringEncode) -> Vec<String> {
    match input {
        StringEncode::Uncoded(vec) => vec,
        StringEncode::Deflate(vec) => {
            let mut decoder = Decoder::new(&vec[..]);
            let mut decoded_data = Vec::new();
            let _result = decoder.read_to_end(&mut decoded_data).unwrap();
            String::from_utf8(decoded_data).unwrap().split('\0').map(|s| s.to_string()).collect()
        }
    }
}
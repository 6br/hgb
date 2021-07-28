use std::io::{self, Read, Result};
use vbyte::{VByteDecoder, VByteEncoded};
//use bitpacking::{BitPacker4x, BitPacker};
use libflate::deflate::{Decoder, Encoder};

pub type VByte = Vec<u8>;
pub type DeltaVByte = Vec<u8>;
pub type Pfor = Vec<u8>;
pub type Deflate = Vec<u8>;

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub enum IntegerEncode {
    Uncoded(Vec<u64>),
    Delta(Vec<u64>),
    VByte(VByte),           // If not sorted
    DeltaVByte(DeltaVByte), // If sorted
    DeltaPfor(Pfor),        // Preserved but not used so far
}

pub enum FloatEncode {
    Uncoded(Vec<f64>),
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub enum StringEncode {
    Uncoded(Vec<String>),
    Deflate(Deflate),
}

/// Returns an encoded array of u64.
pub fn integer_encode(input: &Vec<u64>, sorted: bool) -> Result<IntegerEncode> {
    let mut vec = Vec::with_capacity(input.len());

    if input.is_empty() {
        if sorted {
            return Ok(IntegerEncode::DeltaVByte(vec));
        } else {
            return Ok(IntegerEncode::VByte(vec));
        }
    }

    VByteEncoded::new(input[0]).write_to(&mut vec)?;

    // Delta-encoding
    if sorted {
        if input.len() > 1 {
            for i in 0..input.len() - 1 {
                // vec.push(input[i+1] - input[i]);
                debug_assert!(input[i + 1] >= input[i]);
                VByteEncoded::new(input[i + 1] - input[i]).write_to(&mut vec)?;
            }
        }
        Ok(IntegerEncode::DeltaVByte(vec))
    //vec
    } else {
        if input.len() > 1 {
            for i in 1..input.len() {
                // vec.push(input[i+1] - input[i]);
                VByteEncoded::new(input[i]).write_to(&mut vec)?;
            }
        }
        Ok(IntegerEncode::VByte(vec))
        //vec
    }
}

/// TODO() Remove this: everytime returns vec<u8>
pub fn integer_encode_wrapper(input: &Vec<u64>, sorted: bool) -> Vec<u8> {
    match integer_encode(input, sorted).unwrap() {
        IntegerEncode::DeltaVByte(vec) => vec,
        IntegerEncode::VByte(vec) => vec,
        _ => panic!("Invalid encoding"),
    }
}

/// Returns a decoded array of u64
pub fn integer_decode(input: IntegerEncode) -> Vec<u64> {
    match input {
        IntegerEncode::DeltaVByte(vec) => {
            let decoded = VByteDecoder::new(vec.as_slice()).collect::<Vec<_>>();

            // Delta-decoding
            let mut vec = Vec::with_capacity(decoded.len());
            let mut previous_value = decoded[0];
            vec.push(previous_value);
            if decoded.len() > 1 {
                for i in 1..decoded.len() {
                    vec.push(decoded[i] + previous_value);
                    previous_value = decoded[i];
                }
            }
            vec
        }
        IntegerEncode::VByte(vec) => VByteDecoder::new(vec.as_slice()).collect::<Vec<_>>(),
        IntegerEncode::Uncoded(vec) => vec,
        _ => panic!("Not implemented!"),
    }
}

pub fn string_encode(input: &Vec<String>) -> Vec<u8> {
    let concatenate_string = input.join("\0");
    // println!("str: {}", concatenate_string);
    let mut encoder = Encoder::new(Vec::new());
    io::copy(&mut concatenate_string.as_bytes(), &mut encoder).unwrap();
    let encoded_data = encoder.finish().into_result().unwrap();
    //StringEncode::Deflate(
    encoded_data
}

pub fn string_decode(input: &StringEncode) -> Vec<String> {
    match input {
        StringEncode::Uncoded(vec) => vec.clone(),
        StringEncode::Deflate(vec) => {
            let mut decoder = Decoder::new(&vec[..]);
            let mut decoded_data = Vec::new();
            let _result = decoder.read_to_end(&mut decoded_data).unwrap();
            String::from_utf8(decoded_data)
                .unwrap()
                .split('\0')
                .map(|s| s.to_string())
                .collect()
        }
    }
}

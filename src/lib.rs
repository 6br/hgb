#![feature(rustc_private)]
#[macro_use]

extern crate log;
extern crate libc;

extern crate serde_json;
extern crate bitpacking;
extern crate regex;

#[no_mangle]
pub extern "C" fn hello_rust() -> *const u8 {
    "Hello, world!\0".as_ptr()
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

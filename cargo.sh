#!/bin/bash 

export CACHE_DIR=`pwd`/cargo
export CARGO_HOME=$CACHE_DIR
export RUSTUP_HOME=$CACHE_DIR

curl --proto '=https' --tlsv1.2 https://sh.rustup.rs -sSf | sh -s -- --no-modify-path -y --profile minimal --default-toolchain nightly

mkdir .cargo

target=`rustup show active-toolchain`

echo -e \#\!/bin/sh\\n $@ \"\$@\" > .cargo/cargo_script.sh
chmod +x .cargo/cargo_script.sh
pwd=`pwd`

if [[ ${target} =~ ^[^-]*-(.*)\ .*$ ]]; then
  echo -e "[target.${BASH_REMATCH[1]}]\nlinker = '$pwd/.cargo/cargo_script.sh'" > .cargo/config
fi


`pwd`/cargo/bin/cargo build --release # C-compiler, library path, 


#!/bin/bash
# USAGE: ./image_wrapper.sh READ_MAX <PARAMS>

mkdir -p dnd
READ_MAX=$1
CMD="./target/debug/hgb"
if [ -e "./target/release/hgb" ]; then
  CMD="./target/release/hgb"
fi
shift
ARGS=$@

$CMD -t1 vis -_ 1000000 $ARGS -o dnd/0.png 2>/dev/null

for i in `seq 1 ${READ_MAX}`
do
  $CMD -t1 vis -_ $((i-1)) $ARGS -o dnd/$i.png -* 2>/dev/null
done

$CMD -t1 vis $ARGS -P -A -o dnd/$((READ_MAX+2)).png 2>/dev/null

if [ ! -e "./build.zip" ]; then
  wget https://github.com/6br/interactive-hgb/releases/download/v0.1.0/build.zip 
  unzip build.zip
fi

mv dnd build/
cd build

ret=`python -c 'import sys; print("%i" % (sys.hexversion<0x03000000))'`
if [ $ret -eq 0 ]; then
    #echo "we require python version <3"
    python -m http.server 8000
else
    #echo "python version is <3"
    python -m SimpleHTTPServer 8000
fi
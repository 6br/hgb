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

$CMD -t1 vis -_ 1000000 $ARGS -o dnd/0.png 2>/dev/null &

for i in `seq 1 ${READ_MAX}`
do
  $CMD -t1 vis -_ $((i-1)) $ARGS -o dnd/$i.png -* 2>/dev/null &
done

$CMD -t1 vis $ARGS -P -A -o dnd/$((READ_MAX+2)).png 2>/dev/null &

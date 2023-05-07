#!/bin/bash

my_file=$1

cat "$my_file" | perl -pe "s/[-yYzZ]//g" | sed '/^[[:space:]]*$/d' > "$my_file"_CG.txt


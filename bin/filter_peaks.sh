#!/usr/bin/env bash

input=$1; shift

awk -F "\t" '$7 > 2 && $9 > 2' $input

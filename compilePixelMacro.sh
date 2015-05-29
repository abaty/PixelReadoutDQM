#!/bin/bash

mkdir -p pdfDir
g++ $1 $(root-config --cflags --libs) -Werror -Wall -O2 -o "${1/%.C/}.exe" 
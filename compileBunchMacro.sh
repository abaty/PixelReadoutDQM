#!/bin/bash

g++ $1 -Werror -Wall -O2 -o "${1/%.C/}.exe" 
#!/usr/bin/env bash

make clean
autoreconf -iv --install
# ./configure CC=gcc
./configure --prefix=$HOME
make

make install
make dist

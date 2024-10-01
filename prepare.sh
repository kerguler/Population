#!/usr/bin/env bash

make clean
autoreconf -iv --install
# ./configure CC=gcc
# ./configure --prefix=$HOME
./configure
make

make install
make dist

#!/bin/bash

cd ../gdspy
gcc -O3 -fPIC -shared -I/usr/include/python2.7 -L/usr/lib -lpython2.7 -o boolext.so boolext.c
g++ -O3 -fPIC -shared -I/usr/include/python2.7 -L/usr/lib -lpython2.7 -o clipper.so clipper.cpp
cd -

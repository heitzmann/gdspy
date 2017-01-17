#!/bin/bash

cd gdspy
gcc -O3 -fPIC -shared -I/usr/include/python3.6m -L/usr/lib -lpython3.6m -o boolext.so boolext.c
g++ -O3 -fPIC -shared -I/usr/include/python3.6m -L/usr/lib -lpython3.6m -o clipper.so clipper.cpp
cd -

#!/usr/bin/python

import sys

infile = open(sys.argv[1])
outfile = open(sys.argv[2], 'w')

for line in infile:
    chunk = line.split('\M')
    lines = chunk[0].split('\r')
    for item in lines:
        outfile.write(item + '\n')

infile.close()
outfile.close()

#!/usr/bin/env python3

import sys

filename = sys.argv[1]
queryString = sys.argv[2]

savedText = []

printSwitch = False

with open(filename, 'r') as infile:
    for line in infile:
        if line.startswith('>'):
            if queryString in line:
                printSwitch = True
            else:
                printSwitch = False
        if printSwitch:
            savedText.append(line)

print(''.join(savedText))

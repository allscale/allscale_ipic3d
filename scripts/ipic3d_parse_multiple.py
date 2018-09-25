#!/usr/bin/env python
import math
import sys
import subprocess
import glob

files = glob.glob(sys.argv[1] + "/*.o")

for f in files:
    splitResult = f.split( "_" )
    nodeInfo = splitResult[len(splitResult) - 1]
    command = "grep 'Throughput:' " + f + " | awk -F':' '{print $2}' | awk -F' ' '{print $1}'"

    (res,err) = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()
    res = res.split()

    output = nodeInfo[0:len(nodeInfo)-7] + " "
    for var in res:
        output = output + var + " "
    print output



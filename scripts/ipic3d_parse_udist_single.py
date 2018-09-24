#!/usr/bin/env python
import math
import sys
import subprocess

command = "grep 'Throughput:' " + sys.argv[1] + " | awk -F':' '{print $2}' | awk -F' ' '{print $1}'"

(res,err) = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()
res = res.split()

count = 0
res1 = ""
res2 = ""
res4 = ""
res8 = ""
res16 = ""
res32 = ""
for i in 0, 1, 2, 3:
	for j in 0, 1, 2, 3, 4, 5:
		res1 += res[i*36 + 0 + j] + " "
		res2 += res[i*36 + 6 + j] + " "
		res4 += res[i*36 + 12 + j] + " "
		res8 += res[i*36 + 18 + j] + " "
		res16 += res[i*36 + 24 + j] + " "
		res32 += res[i*36 + 30 + j] + " "

print res1
print res2
print res4
print res8
print res16
print res32

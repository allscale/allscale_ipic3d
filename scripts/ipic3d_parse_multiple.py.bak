#!/usr/bin/env python
import math
import sys
import subprocess

command = "grep 'Throughput:' " + sys.argv[1] + " | awk -F':' '{print $2}' | awk -F' ' '{print $1}'"

(res,err) = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()
res = res.split()

count = 0
sum = 0
cores = int(sys.argv[2])
output = ""
for var in res:
	sum = sum + float(var)
	count += 1
	if count % cores == 0:
		#print sum / cores
		#'{:.2e}'.format(sum / cores)
		output += '{0:2.5e}'.format(sum / cores) + "  " 
		count = 0
		sum = 0

print output



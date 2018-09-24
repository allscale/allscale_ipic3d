
TMP1=$(grep "Throughput:" allscale_ipic3d_512K_1node.o | awk -F':' '{print $2}' | awk -F' ' '{print $1}')
echo $TMP1
echo $TMP1[0]

TMP=$(grep "Throughput:" allscale_ipic3d_512K_2nodes.o | awk -F':' '{print $2}' | awk -F' ' '{print $1}')
echo $TMP
#sum=0
#count=1
#for var in "${TMP[@]}"
#do
#	sum=$((sum+var))
#	((count++))
#done
#echo $sum
#echo $count

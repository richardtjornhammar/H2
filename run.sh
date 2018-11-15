#!/bin/bash
x=0;y=0;z=0;
echo " "> outpTE

for k in `seq 0.15 0.01 4.00`;

do

# SIMULATE ATT DIFFERENT SEPARATIONS

#echo "0 3 2 1 0" > inp
#echo "$k 2.0925 1.24" >> inp

echo "0 3 1 1 0" > inp
echo "$k 2.0925 1.24" >> inp

./quantum 2> temp
cat temp|grep "Total energy">> outpTE
cat temp|grep "Electronic energy"> outpEE

done

#1.4632 2.0925 1.24

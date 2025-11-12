#!/bin/bash

for i in *.bw
do
 /home/martinalavanya/Documents/ETH/mls_semesterprojekt/kentutils/bigWigToBedGraph $i ${i/%bw/bedGraph}
done


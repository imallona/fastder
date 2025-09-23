#!/bin/bash

for i in *.bw
do
 ../../utilities/bigWigToBedGraph $i ${i/%bw/bedGraph}
done


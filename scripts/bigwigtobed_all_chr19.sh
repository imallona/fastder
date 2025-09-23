#!/bin/bash

for i in *.bw
do
 ../../utilities/bigWigToBedGraph -chrom=chr19 $i ${i/%.bw/_chr19.bedGraph}
done


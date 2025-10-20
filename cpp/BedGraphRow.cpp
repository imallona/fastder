//
// Created by marti on 08/10/2025.
//
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "BedGraphRow.h"

// print a BedGraphRow
void BedGraphRow::print() const {
    std::cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << total_reads <<  "\t" << length << std::endl;
}


//normalize reads to CPM for better comparability between libraries
void BedGraphRow::normalize(const double& library_size)
{
    this->coverage = (this->coverage / library_size) * 1e6;
    //coverage is NOT cumulative across the bin but rather per base pair coverage within the bin
}
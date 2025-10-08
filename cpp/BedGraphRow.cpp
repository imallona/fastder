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

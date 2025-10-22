//
// Created by marti on 08/10/2025.
//
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <BedGraphRow.h>
#include <cstdint>


// constructor with arguments
BedGraphRow::BedGraphRow(std::string _chrom, uint64_t _start, uint64_t _end, double _coverage) : chrom(_chrom),
    start(_start), end(_end), coverage(_coverage), total_reads(0), length(_end - _start)
{
}

// print a BedGraphRow
void BedGraphRow::print() const {
    std::cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t";
    if (total_reads > 0)
        std::cout << total_reads <<  "\t";
    std::cout << length << std::endl;
}


//normalize reads to CPM for better comparability between libraries
void BedGraphRow::normalize(const uint64_t library_size)
{
    this->coverage = (this->coverage / library_size) * 1e6;
    //coverage is NOT cumulative across the bin but rather per base pair coverage within the bin
}
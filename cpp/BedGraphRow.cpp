//
// Created by marti on 08/10/2025.
//
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <BedGraphRow.h>
#include <cstdint>


// constructor with arguments. Coverage rows are unstranded by default since the
// monorail-style BigWig the workflow emits is not strand-resolved.
BedGraphRow::BedGraphRow(std::string _chrom, uint64_t _start, uint64_t _end, double _coverage) : chrom(_chrom),
    start(_start), end(_end), coverage(_coverage), total_reads(0), length(_end - _start), strand('.')
{
}

// constructor for stranded BigWig input. _strand is '+' or '-'; '.' means
// unstranded (the BigWig contained reads from both strands).
BedGraphRow::BedGraphRow(std::string _chrom, uint64_t _start, uint64_t _end, double _coverage, char _strand) :
    chrom(_chrom), start(_start), end(_end), coverage(_coverage), total_reads(0),
    length(_end - _start), strand(_strand)
{
}

// prints a BedGraphRow
void BedGraphRow::print() const {
    std::cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t";
    if (total_reads > 0)
        std::cout << total_reads <<  "\t";
    std::cout << length << std::endl;
}


// normalizes reads to CPM for better comparability
void BedGraphRow::normalize(const uint64_t library_size)
{
    if (library_size == 0)
    {
        std::cerr << "[ERROR] Library size is zero." << std::endl;
        return;
    }
    this->coverage = (this->coverage / library_size) * 1e6;
    //coverage is NOT cumulative across the bin but rather per base pair coverage within the bin
}
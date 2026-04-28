//
// Created by marti on 08/10/2025.
//
#pragma once

#ifndef FASTDER_BEDGRAPHROW_H
#define FASTDER_BEDGRAPHROW_H

#include <string>
#include <cstdint>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <BedGraphRow.h>
#include <cstdint>

#endif //FASTDER_BEDGRAPHROW_H

class BedGraphRow
{
public:
    std::string chrom;
    uint64_t start;
    uint64_t end;
    double coverage; // normalized coverage by CPM
    unsigned int total_reads; // number of reads spanning across the bin, total_reads = length * coverage
    unsigned int length;
    // strand is '.' for unstranded coverage (the default), '+' for plus-strand
    // BigWigs and '-' for minus-strand BigWigs. Keeping the field on every row
    // means stranded data can flow through the pipeline without a parallel
    // data structure. Downstream code currently ignores it; future work in
    // find_ERs / stitch_up can filter on it.
    char strand = '.';
    // add optional values for average coverage, DER identifier


    BedGraphRow() = default;
    BedGraphRow(std::string chrom, uint64_t start, uint64_t end, double coverage);
    BedGraphRow(std::string chrom, uint64_t start, uint64_t end, double coverage, char strand);
    void print() const;
    void normalize(const uint64_t library_size);

};
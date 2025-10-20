//
// Created by marti on 08/10/2025.
//
#pragma once

#ifndef FASTDER_BEDGRAPHROW_H
#define FASTDER_BEDGRAPHROW_H

#endif //FASTDER_BEDGRAPHROW_H

class BedGraphRow
{
public:
    std::string chrom;
    unsigned int start;
    unsigned int end;
    double coverage; // normalized coverage by CPM
    unsigned int total_reads; // number of reads spanning across the bin, total_reads = length * coverage
    unsigned int length;
    // add optional values for average coverage, DER identifier

    //print BedGraphRow
    BedGraphRow() = default;
    void print() const;
    void normalize(const double& library_size);

};
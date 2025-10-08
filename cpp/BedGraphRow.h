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
    double coverage;
    //double avg = 0;
    unsigned int total_reads;
    unsigned int length;
    // add optional values for average coverage, DER identifier

    //print BedGraphRow
    BedGraphRow() = default;
    void print() const;

};
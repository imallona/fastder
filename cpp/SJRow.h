//
// Created by marti on 08/10/2025.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#ifndef FASTDER_SPLICE_JUNCTION_H
#define FASTDER_SPLICE_JUNCTION_H
#include <cstdint>

class SJRow
{
public:
    std::string chrom;
    uint64_t start;
    uint64_t end;
    unsigned int length;
    bool strand; // 1 = +, 0 = -
    bool annotated; // 0 or 1

    // The recount3 RR file carries four extra per-junction columns
    // (left_motif, right_motif, left_annotated, right_annotated) that fastder
    // reads but never accesses after parse. On full hg38 with about 9.5M
    // junctions those strings cost roughly 2 GB of resident heap. operator>>
    // below parses and discards them so the on-disk format stays compatible
    // and the storage cost goes away.

    // constructor
    SJRow() = default;
    SJRow(std::string _chrom, uint64_t _start, uint64_t _end, int _length, char _strand, bool _annotated);

    // overload input operator for SJRow
    friend std::istream& operator>>(std::istream &is, SJRow &row) {
        char strand_;
        std::string discard;
        if (!(is >> row.chrom >> row.start >> row.end >> row.length >> strand_ >> row.annotated
                 >> discard >> discard >> discard >> discard)) {
            return is;
        }
        row.strand = (strand_ == '+'); // 1 if +
        return is;
    }


    // overload output operator for SJRow (debug only, not part of fastder's outputs)
    friend std::ostream& operator<< (std::ostream& os, const SJRow& row)
    {
        return os << row.chrom << "\t" << row.start << "\t" << row.end << "\t" << row.length
                  << "\t" << (row.strand ? '+' : '-') << "\t" << row.annotated;
    }

};




#endif //FASTDER_SPLICE_JUNCTION_H

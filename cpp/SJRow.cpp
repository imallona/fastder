//
// Created by marti on 08/10/2025.
//

#include "SJRow.h"

SJRow::SJRow(std::string _chrom, uint64_t _start, uint64_t _end, int _length, char _strand, bool _annotated) {
    chrom = _chrom;
    start = _start;
    end = _end;
    length = _length;
    strand = (_strand == '+');
    annotated = _annotated;
};

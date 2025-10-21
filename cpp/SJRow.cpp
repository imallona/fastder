//
// Created by marti on 08/10/2025.
//

#include "SJRow.h"

SJRow::SJRow(std::string _chrom, uint64_t _start, uint64_t _end, int _length, char _strand, bool _annotated,
             std::string _left_motif, std::string _right_motif, std::string _left_annotated, std::string _right_annotated) {
    chrom = _chrom;
    start = _start;
    end = _end;
    length = _length;
    strand = (_strand == '+');
    annotated = _annotated;
    left_motif = _left_motif;
    right_motif = _right_motif;
    left_annotated = _left_annotated;
    right_annotated = _right_annotated;
};


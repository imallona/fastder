//
// Created by marti on 08/10/2025.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#ifndef FASTDER_SPLICE_JUNCTION_H
#define FASTDER_SPLICE_JUNCTION_H


class SJRow
{
public:
    std::string chrom;
    int start;
    int end;
    int length;
    bool strand; // 1 = +, 0 = -
    bool annotated; // 0 or 1
    std::string left_motif;
    std::string right_motif;
    std::string  left_annotated;
    std::string  right_annotated;

    // constructor
    SJRow() = default;
    SJRow(std::string _chrom, int _start, int _end, int _length, char _strand, bool _annotated,
        std::string _left_motif, std::string _right_motif, std::string _left_annotated, std::string _right_annotated);

    // overload input operator for SJRow
    friend std::istream& operator>>(std::istream &is, SJRow &row) {
        char strand_;
        is >> row.chrom >> row.start >> row.end >> row.length >> strand_ >> row.annotated >> row.left_motif
        >> row.right_motif >> row.left_annotated >> row.right_annotated;
        row.strand = (strand_ == '+'); // 1 if +
        return is;
    }


    // overload output operator for SJRow
    friend std::ostream& operator<< (std::ostream& os, const SJRow& row)
    {
        return os << row.chrom << "\t" << row.start << "\t" << row.end << "\t" << row.length << "\t" << row.strand << "\t" <<
            row.annotated << "\t" << row.left_motif << "\t" << row.right_motif << "\t" << row.left_annotated << "\t" << row.right_annotated;

    }

};




#endif //FASTDER_SPLICE_JUNCTION_H
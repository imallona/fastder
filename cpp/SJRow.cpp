//
// Created by marti on 08/10/2025.
//

#include "SJRow.h"

SJRow::SJRow(std::string _chrom, int _start, int _end, int _length, char _strand, bool _annotated,
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
// overload input operator for SJRow
std::istream &SJRow::operator>>(std::istream &is) {
    char annotated;
    is >> this.chrom >> row.start >> row.end >> row.length >> row.strand >> annotated >> row.left_motif
    >> row.right_motif >> row.left_annotated >> row.right_annotated;
    row.annotated = (annotated == '+'); // 1 if +
    return is;
}


// overload output operator for SJRow
std::ostream& SJRow::operator<< (std::ostream& os)
{
    return os << chrom << "\t" << start << "\t" << end << "\t" << length << "\t" << strand << "\t" <<
        annotated << "\t" << left_motif << "\t" << right_motif << "\t" << left_annotated << "\t" << right_annotated;

}
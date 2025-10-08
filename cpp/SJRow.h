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
    SJRow(std::string _chrom, int _start, int _end, int _length, char _strand, bool _annotated,
        std::string _left_motif, std::string _right_motif, std::string _left_annotated, std::string _right_annotated);

    std::ostream& operator<< (std::ostream& os);


};



#endif //FASTDER_SPLICE_JUNCTION_H
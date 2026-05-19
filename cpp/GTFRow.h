//
// Created by martinalavanya on 28.10.25.
//

#ifndef MLS_GTF_H
#define MLS_GTF_H
#include "StitchedER.h"


class GTFRow
{

    public:
    GTFRow(const StitchedER& region, std::string ftr, unsigned int id);
    void change_feature(std::string ftr, unsigned int id, unsigned int exon_nr);

    // MEMBER VARIABLES
    std::string seqname; //name of the chromosome without the "chr" prefix
    std::string source = "fastder";
    std::string feature;
    uint64_t start;
    uint64_t end;
    double score; // the avg coverage of the stitched ER
    std::string strand = ".";
    std::string frame = ".";
    std::string attribute;

    // overload output operator for GTFRow
    friend std::ostream& operator<< (std::ostream& os, const GTFRow& row)
    {
        // Internal coordinates are 0-based half-open (BedGraph / BigWig
        // convention). GTF is 1-based inclusive, so the start is shifted by
        // one. The end needs no shift: a 0-based half-open end equals the
        // 1-based inclusive end of the same interval.
        return os << row.seqname << "\t" << row.source << "\t" << row.feature << "\t" << (row.start + 1) << "\t" << row.end << "\t" <<
            row.score << "\t" << row.strand << "\t" << row.frame << "\t" << row.attribute;

    }
};



#endif //MLS_GTF_H
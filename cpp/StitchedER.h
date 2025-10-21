//
// Created by martinalavanya on 20.10.25.
//

#ifndef MLS_STITCHEDER_H
#define MLS_STITCHEDER_H

#endif //MLS_STITCHEDER_H

#include <vector>
#include <BedGraphRow.h>
#include <cstdint>

class StitchedER
{
public:
    // MEMBER FUNCTIONS
    StitchedER() = default;

    StitchedER(const BedGraphRow& expressed_region, unsigned int er_id);
    void append(unsigned int er_id,unsigned int length, double coverage);
    double get_avg_coverage();

    //bool is_similar(double val1, double val2);


    // MEMBER VARIABLES

    std::vector<unsigned int> er_ids; //all expressed regions in a stitched_ER, er_id corresponds to index of averager.expressed_regions
    // example: stitched_ER consists of er_ids 45, 46, 47, 49 == vector indices of expressed_regions
    double across_er_coverage; // avg (weighted) coverage of all exons that are part of the stitched ER so far
    std::vector<std::pair<unsigned int, double>> all_coverages; // stores a pair of er length (= weight) + normalized average coverage of the er
    unsigned int total_length;
    uint64_t start; //TODO convert to uint64_t?
    uint64_t end;


    // overload output operator for SJRow
    friend std::ostream& operator<< (std::ostream& os, const StitchedER& stitched_er)
    {

        for (unsigned int i = 0; i < stitched_er.er_ids.size(); i++)
        {
            os << i << "\t" << "(" <<  stitched_er.all_coverages[i].first <<"," << stitched_er.all_coverages[i].second << ")" << std::endl;
        }
        return os << stitched_er.across_er_coverage << "\t" << stitched_er.start << "\t" << stitched_er.end << "\t" << stitched_er.total_length << std::endl;

    }
};
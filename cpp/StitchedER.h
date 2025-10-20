//
// Created by martinalavanya on 20.10.25.
//

#ifndef MLS_STITCHEDER_H
#define MLS_STITCHEDER_H

#endif //MLS_STITCHEDER_H

#include <vector>

class StitchedER
{
    // MEMBER FUNCTIONS
    double get_avg_coverage();
    bool is_similar(double val1, double val2);


    // MEMBER VARIABLES

    std::vector<unsigned int> er_ids; //all expressed regions in a stitched_ER, er_id corresponds to index of averager.expressed_regions
    // example: stitched_ER consists of er_ids 45, 46, 47, 49 == vector indices of expressed_regions
    double across_er_coverage; // avg (weighted) coverage of all exons that are part of the stitched ER so far
    std::vector<std::pair<int, double>> all_coverages; // stores a pair of er length (= weight) + normalized average coverage of the er
    unsigned int total_reads = 0;
    unsigned int length = 0;
};
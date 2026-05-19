//
// Created by martinalavanya on 20.10.25.
//

#include "StitchedER.h"


// initialize a new StitchedER from one ER
StitchedER::StitchedER(const BedGraphRow& expressed_region, int er_id)
{
    chrom = expressed_region.chrom;
    start = expressed_region.start;
    end = expressed_region.end;
    er_ids = {er_id};
    across_er_coverage = expressed_region.coverage;
    all_coverages = {std::make_pair(expressed_region.length, expressed_region.coverage)};
    er_bounds = {std::make_pair(expressed_region.start, expressed_region.end)};
    total_length = expressed_region.length;
}

// later change to add weight / memory
double StitchedER::get_avg_coverage() const
{
    double sum = 0;
    unsigned int full_length = 0;
    // all_coverages: <length, coverage>
    for (auto& er : this->all_coverages){
        // actually should be 0.0, but using a slightly larger value to avoid exact comparison of doubles
        if (er.second >= 0.000001)
        {
            sum += (er.first * er.second); // weighted sum of avg coverages across ERs
            full_length += er.first;
        }
    }
    if (full_length == 0)
    {
        std::cerr << "[ERROR] All stitched ERs have a combined length of 0. This shouldn't happen..." << std::endl;
        return sum;
    }
    return sum / full_length;
}

// note that er_id is 0-indexed. lo/hi is the ER's coverage extent, or the
// splice junction's [donor, acceptor] when er_id is -1 (a spliced region).
void StitchedER::append(int er_id, unsigned int length, double coverage,
                        uint64_t lo, uint64_t hi)
{
    er_ids.emplace_back(er_id); // -1 for spliced regions
    all_coverages.push_back({std::make_pair(length, coverage)});
    er_bounds.push_back(std::make_pair(lo, hi));
    // only update avg coverage if an exon was added
    if (er_id != -1)
    {
        across_er_coverage = this->get_avg_coverage();
        total_length += length; // only count er length
    }
    end += length; // count er_length and sj_length to get the correct end position


}




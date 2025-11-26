//
// Created by martinalavanya on 20.10.25.
//

#include "StitchedER.h"


// initialize a new StitchedER with one ER in it
StitchedER::StitchedER(const BedGraphRow& expressed_region, int er_id)
{
    chrom = expressed_region.chrom;
    start = expressed_region.start;
    end = expressed_region.end;
    er_ids = {er_id};
    across_er_coverage = expressed_region.coverage;
    all_coverages = {std::make_pair(expressed_region.length, expressed_region.coverage)};
    total_length = expressed_region.length;
}

// later change to add weight / memory
double StitchedER::get_avg_coverage(){
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
    return sum / full_length;
}

// note that er_id is 0-indexed
void StitchedER::append(int er_id, unsigned int length, double coverage)
{

    er_ids.push_back(er_id); // -1 for spliced regions
    all_coverages.push_back({std::make_pair(length, coverage)});
    // only update avg coverage if an exon was added
    if (er_id != -1)
    {
        across_er_coverage = this->get_avg_coverage();
        total_length += length; // only count er length
    }

    //std::cout << "Adding length = " << length << " to position = " << end << std::endl ;
    end += length; // count er_length and sj_length to get the correct end position


}




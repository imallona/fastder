//
// Created by martinalavanya on 20.10.25.
//

#include "StitchedER.h"


// initialize a new StitchedER with one ER in it
StitchedER::StitchedER(const BedGraphRow& expressed_region, unsigned int er_id)
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
    for (auto& er : this->all_coverages){
        sum += (er.first * er.second); // weighted sum of avg coverages across ERs
        full_length += er.first;

    }
    return sum / full_length;
}

// note that er_id is 0-indexed
void StitchedER::append(unsigned int er_id, unsigned int sj_length, unsigned int er_length, double coverage)
{

    er_ids.push_back(er_id);
    all_coverages.push_back({std::make_pair(er_length, coverage)});
    across_er_coverage = this->get_avg_coverage();
    total_length += er_length; // only count er length
    end = end + sj_length + er_length; // count er_length and sj_length to get the correct end position

}




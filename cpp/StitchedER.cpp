//
// Created by martinalavanya on 20.10.25.
//

#include "StitchedER.h"


// initialize a new StitchedER with one ER in it
//TODO maybe just pass the BedGraphRow object instead
StitchedER::StitchedER(uint64_t start_, uint64_t end_, int er_id_, unsigned int length_, double coverage_)
{
    start = start_;
    end = end_;
    er_ids = {er_id_};
    across_er_coverage = coverage_;
    all_coverages = {std::make_pair(length_, coverage_)};
    total_length = length_;
}

void StitchedER::append(unsigned int er_id,unsigned int length, double coverage)
{
    er_ids.push_back(er_id);
    all_coverages.push_back({std::make_pair(length, coverage)});
    across_er_coverage = get_avg_coverage();
    total_length += length;

}
// later change to add weight / memory
double StitchedER::get_avg_coverage(){
    double sum = 0;
    unsigned int full_length = 0;
    for (auto& er : this->all_coverages){
        sum += er.first * er.second; // weighted sum of avg coverages across ERs
        full_length += er.first;

    }
    return sum / total_length;
}



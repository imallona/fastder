//
// Created by martinalavanya on 20.10.25.
//

#include "Integrator.h"
#include "BedGraphRow.h"
#include "StitchedER.h"
#include "SJRow.h"

void Integrator::stitch_up(const std::vector<BedGraphRow>& expressed_regions, const std::unordered_map<unsigned int, unsigned int>& mm_sj_counts, const std::vector<SJRow>& rr)
{

    StitchedER er1; // define the first StitchedER, currently consisting of 1 ER
    stitched_ERs.push_back();
        // iterate over expressed regions
    for (auto& region : expressed_regions)
    {
      auto current_sj = mm_sj_counts.begin()->first; // stores first sj id


    }

    }
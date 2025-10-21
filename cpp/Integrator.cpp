//
// Created by martinalavanya on 20.10.25.
//

#include <Integrator.h>

// function that calculates relative match with a tolerance of +/- 5%
bool Integrator::within_threshold(double val1, double val2){
    double tolerance_bottom = val1 * 0.95;
    double tolerance_top = val1 * 1.05;
    return val2 >= tolerance_bottom && val2 <= tolerance_top;
}

bool Integrator::within_threshold(uint64_t val1, uint64_t val2){
    double tolerance_bottom = val1 * 0.95;
    double tolerance_top = val1 * 1.05;
    return val2 >= tolerance_bottom && val2 <= tolerance_top;
}

// function that calculates relative match with a tolerance of +/- 5%
bool Integrator::is_similar(const StitchedER& most_recent_er, const BedGraphRow& expressed_region, const SJRow& current_sj){

    return (within_threshold(most_recent_er.end, current_sj.start)
       && within_threshold(expressed_region.start, current_sj.end)
       && within_threshold(most_recent_er.across_er_coverage, expressed_region.coverage)); //TODO maybe compare with across_er_coverage instead
}

void Integrator::stitch_up(const std::vector<BedGraphRow>& expressed_regions, const std::unordered_map<unsigned int, unsigned int>& mm_sj_counts, const std::vector<SJRow>& rr_all_sj)
{

    StitchedER er1 = StitchedER(expressed_regions[0].start, expressed_regions[0].end, 0, expressed_regions[0].length, expressed_regions[0].coverage); // define the first StitchedER, currently consisting of 1 ER

    stitched_ERs.push_back(er1); //TODO MAYBE NOT YET, ONLY APPEND WHEN IT'S FINISHED
    auto current_sj = mm_sj_counts.begin(); //iterator over the vector of SJ ids
        // iterate over expressed regions
    for (unsigned int i = 0; i < expressed_regions.size(); ++i)
    {
        auto expressed_region = expressed_regions[i];
        StitchedER& most_recent_er = stitched_ERs.back(); // this is one expressed region right now

        // get rr_all_sj, which is a vector of SJRows
        if (is_similar(most_recent_er, expressed_region, rr_all_sj[current_sj->first]))
        {
            most_recent_er.append(i, expressed_region.length, expressed_region.coverage);
            // move to next SJ
            ++current_sj;
        }

        else if ()


    }

    }
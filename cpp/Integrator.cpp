//
// Created by martinalavanya on 20.10.25.
//

#include <Integrator.h>
// constructor
Integrator::Integrator()
{
    stitched_ERs = {}; // empty vector with elements of type StitchedER
}



// function that calculates relative match with a tolerance of +/- 5%
bool Integrator::within_threshold(double val1, double val2){
    double tolerance_bottom = val1 * (1 - tolerance);
    double tolerance_top = val1 * (1 + tolerance);
    return val2 >= tolerance_bottom && val2 <= tolerance_top;
}

bool Integrator::within_threshold(uint64_t val1, uint64_t val2){
    double tolerance_bottom = val1 * (1 - tolerance);
    double tolerance_top = val1 * (1 + tolerance);
    return val2 >= tolerance_bottom && val2 <= tolerance_top;
}

// function that calculates relative match with a tolerance of +/- 5%
bool Integrator::is_similar(const StitchedER& most_recent_er, const BedGraphRow& expressed_region, const SJRow& current_sj){

    return (within_threshold(most_recent_er.end, current_sj.start)
       && within_threshold(expressed_region.start, current_sj.end)
       && within_threshold(most_recent_er.across_er_coverage, expressed_region.coverage)); //TODO maybe compare with across_er_coverage instead
}

// function that calculates relative match with a tolerance of +/- 5%
bool Integrator::sj_too_far_back(const uint64_t most_recent_er_end, const uint64_t sj_start){

    return most_recent_er_end > sj_start
    && !within_threshold(most_recent_er_end, sj_start);
}

void Integrator::stitch_up(const std::vector<BedGraphRow>& expressed_regions, const std::map<uint64_t, unsigned int>& mm_sj_counts, const std::vector<SJRow>& rr_all_sj)
{

    StitchedER er1 = StitchedER(expressed_regions[0], 0); // define the first StitchedER, currently consisting of 1 ER
    stitched_ERs.push_back(er1); //TODO MAYBE NOT YET, ONLY APPEND WHEN IT'S FINISHED
    auto current_sj = mm_sj_counts.begin(); //iterator over the vector of SJ ids
    std::cout << stitched_ERs.front() << std::endl;
    int max_stitched_ers = 0;
    // iterate over expressed regions
    for (unsigned int i = 0; i < expressed_regions.size(); ++i)
    {
        //only compare if we aren't at the last SJ yet
        if (current_sj != mm_sj_counts.end()){

            const auto& expressed_region = expressed_regions[i];
            StitchedER& most_recent_er = stitched_ERs.back(); // this is one expressed region right now


            while (current_sj != mm_sj_counts.end() && sj_too_far_back(most_recent_er.end, rr_all_sj[current_sj->first].start))
            {
                ++current_sj;
            }
            // get rr_all_sj, which is a vector of SJRows
            if (is_similar(most_recent_er, expressed_region, rr_all_sj[current_sj->first]))
            {
                //expressed_region.print();
                most_recent_er.append(i, expressed_region.length, expressed_region.coverage);
                std::cout << "STITCHED region: current er_id = " << i << std::endl;
                std::cout << stitched_ERs.back() << std::endl;
                // move to next SJ
                ++current_sj;

                // find maximum number of ERs that were stitched together
                if (max_stitched_ers < stitched_ERs.back().er_ids.size())
                {
                    max_stitched_ers = stitched_ERs.back().er_ids.size();
                }
            }


            // current ER doesn't belong to any existing ERs --> start a new ER
            else
            {
                stitched_ERs.push_back(StitchedER(expressed_region, i));
            }


        }
    }
    std::cout << "max_stitched_ers =" << max_stitched_ers<< std::endl;

    }
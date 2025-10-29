//
// Created by martinalavanya on 20.10.25.
//

#ifndef MLS_INTEGRATOR_H
#define MLS_INTEGRATOR_H
#include <unordered_map>

#endif //MLS_INTEGRATOR_H


#include <SJRow.h>
#include <BedGraphRow.h>
#include <chrono>
#include <GTFRow.h>
#include <map>

class Integrator
{
    public:
    Integrator();

    void stitch_up(const std::vector<BedGraphRow>& expressed_regions, const std::map<uint64_t, unsigned int>& mm_sj_counts, const std::vector<SJRow>& rr_all_sj);
    bool within_threshold(double val1, double val2);
    bool within_threshold(uint64_t val1, uint64_t val2);
    bool is_similar(const StitchedER& most_recent_er, const BedGraphRow& expressed_region, const SJRow& current_sj);
    bool sj_too_far_back(uint64_t most_recent_er_end, uint64_t sj_start);

    void write_to_gtf(const std::string& output_path);

    // MEMBERS
    std::vector<StitchedER> stitched_ERs;

    double coverage_tolerance = 0.1;
    int position_tolerance = 5;


};
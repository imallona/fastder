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
#include <unordered_set>

class Integrator
{
    public:
    Integrator(double coverage_tolerance_, int position_tolerance_);

    void stitch_up(std::unordered_map<std::string, std::vector<BedGraphRow>>& expressed_regions, const std::unordered_map<std::string, std::vector<uint32_t>>& mm_chrom_sj, const std::vector<SJRow>& rr_all_sj);
    // Per-strand pass invoked by stitch_up. Builds chains of ERs connected
    // by SJs of the requested strand and appends chains of length 2 or
    // more to output, tagged with strand. Indices of ERs that ended up
    // in any emitted chain (including ERs claimed by an earlier-running
    // strand pass) are tracked in consumed_indices, so the second pass
    // skips ERs the first one already used and stitch_up can emit the
    // remaining ERs as unstranded standalones afterwards.
    void stitch_one_strand(const std::string& chrom,
                           char strand,
                           const std::vector<uint32_t>& strand_sjs,
                           const std::vector<BedGraphRow>& ers,
                           const std::vector<SJRow>& rr_all_sj,
                           std::unordered_set<int>& consumed_indices,
                           std::vector<StitchedER>& output);
    bool within_threshold(double val1, double val2) const;
    bool within_threshold(uint64_t val1, uint64_t val2) const;
    bool is_similar(const StitchedER& most_recent_er, const BedGraphRow& expressed_region, const SJRow& current_sj);
    bool sj_too_far_back(uint64_t most_recent_er_end, uint64_t sj_start);

    void write_to_gtf(const std::string& output_path);

    // MEMBERS
    std::vector<StitchedER> stitched_ERs;

    double coverage_tolerance = 0.1;
    int position_tolerance = 5;


};
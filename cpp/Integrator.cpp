//
// Created by martinalavanya on 20.10.25.
//

#include <cassert>
#include <algorithm>
#include <Integrator.h>

// constructor
Integrator::Integrator(double coverage_tolerance_, int position_tolerance_)
{
    stitched_ERs = {}; // empty vector with elements of type StitchedER
    coverage_tolerance = coverage_tolerance_;
    position_tolerance = position_tolerance_;
}

// function that calculates relative match with a tolerance of +/- n%
bool Integrator::within_threshold(double val1, double val2) const
{
    double tolerance_bottom = val1 * (1 - coverage_tolerance);
    double tolerance_top = val1 * (1 + coverage_tolerance);
    return val2 >= tolerance_bottom && val2 <= tolerance_top;
}

// same function for 64-bit integers
bool Integrator::within_threshold(uint64_t pos_1, uint64_t pos_2) const
{
    return pos_1 >= pos_2 - position_tolerance && pos_1 <= pos_2 + position_tolerance;
}

// check if SJ and ERs match (coordinate and coverage check)
bool Integrator::is_similar(const StitchedER& most_recent_er, const BedGraphRow& expressed_region, const SJRow& current_sj){

    return (within_threshold(most_recent_er.end, current_sj.start)
       && within_threshold(expressed_region.start, current_sj.end)
       && within_threshold(most_recent_er.across_er_coverage, expressed_region.coverage)); //TODO perhaps compare with across_er_coverage instead
}

// function that calculates relative match with a tolerance of +/- 5%
bool Integrator::sj_too_far_back(const uint64_t most_recent_er_end, const uint64_t sj_start){

    return most_recent_er_end > sj_start
    && !within_threshold(most_recent_er_end, sj_start);
}

// Per-strand chain builder. Walks the chromosome's ERs left to right, trying
// to extend a candidate chain via SJs of the requested strand. Emits chains
// of two or more ERs to output, tagged with strand. ERs already claimed by
// an earlier strand pass (present in consumed_indices on entry) are skipped:
// they cannot seed a new candidate, they cannot extend a candidate, and
// hitting one closes the current candidate. Indices of ERs that wind up in
// an emitted chain are added to consumed_indices.
void Integrator::stitch_one_strand(const std::string& chrom,
                                   char strand,
                                   const std::vector<uint32_t>& strand_sjs,
                                   const std::vector<BedGraphRow>& ers,
                                   const std::vector<SJRow>& rr_all_sj,
                                   std::unordered_set<int>& consumed_indices,
                                   std::vector<StitchedER>& output)
{
    if (ers.empty() || strand_sjs.empty()) return;

    auto sj_it = strand_sjs.begin();

    // find first non-consumed ER for the seed
    int seed_i = 0;
    while (seed_i < static_cast<int>(ers.size()) && consumed_indices.count(seed_i)) ++seed_i;
    if (seed_i >= static_cast<int>(ers.size())) return;

    StitchedER candidate(ers[seed_i], seed_i);
    candidate.strand = strand;
    bool have_candidate = true;

    auto close_candidate = [&]()
    {
        if (candidate.er_ids.size() > 1)
        {
            output.emplace_back(candidate);
            for (int id : candidate.er_ids)
            {
                if (id >= 0) consumed_indices.insert(id);
            }
        }
        have_candidate = false;
    };

    int max_chain_len = 0;

    for (int i = seed_i + 1; i < static_cast<int>(ers.size()); ++i)
    {
        if (consumed_indices.count(i))
        {
            // ER was claimed by an earlier strand pass; close any in-progress chain
            close_candidate();
            continue;
        }

        const auto& expressed_region = ers[i];
        if (!have_candidate)
        {
            candidate = StitchedER(expressed_region, i);
            candidate.strand = strand;
            have_candidate = true;
            continue;
        }

        if (sj_it == strand_sjs.end())
        {
            close_candidate();
            candidate = StitchedER(expressed_region, i);
            candidate.strand = strand;
            have_candidate = true;
            continue;
        }

        // skip past SJs that lie too far behind the chain's end
        while (sj_it != strand_sjs.end()
            && (candidate.end > rr_all_sj[*sj_it - 1].start
                && !within_threshold(candidate.end, rr_all_sj[*sj_it - 1].start))
            && rr_all_sj[*sj_it - 1].chrom == chrom)
        {
            ++sj_it;
        }
        const auto sj_to_check = (sj_it == strand_sjs.end()) ? std::prev(strand_sjs.end()) : sj_it;

        if (is_similar(candidate, expressed_region, rr_all_sj[*sj_to_check - 1]))
        {
            uint64_t sj_length = expressed_region.start - ers[candidate.er_ids.back()].end;
            candidate.append(-1, sj_length, 0.0);
            candidate.append(i, expressed_region.length, expressed_region.coverage);
            int er_count = 0;
            for (int id : candidate.er_ids) if (id >= 0) ++er_count;
            if (er_count > max_chain_len) max_chain_len = er_count;
            if (sj_it != strand_sjs.end()) ++sj_it;
        }
        else
        {
            close_candidate();
            candidate = StitchedER(expressed_region, i);
            candidate.strand = strand;
            have_candidate = true;
        }
    }

    if (have_candidate) close_candidate();

    if (max_chain_len > 0)
    {
        std::cout << "[INFO] Longest stitched ER in " << chrom << " (" << strand << ") contains "
                  << max_chain_len << " ERs" << std::endl;
    }
}


void Integrator::stitch_up(std::unordered_map<std::string, std::vector<BedGraphRow>>& expressed_regions, const std::unordered_map<std::string, std::vector<uint32_t>>& mm_chrom_sj, const std::vector<SJRow>& rr_all_sj)
{
    // Strand-aware stitching. For each chromosome with at least one ER:
    //   1. Bucket SJs by strand (SJRow.strand: true -> '+', false -> '-').
    //      A chromosome with no SJs at all skips both strand passes and
    //      every ER is emitted as a single-ER StitchedER with strand '.'.
    //   2. Run stitch_one_strand for each non-empty bucket. The shared
    //      consumed_indices set ensures an ER can be in at most one chain.
    //   3. Emit any ER not pulled into a chain as a single-ER StitchedER
    //      with strand '.'.
    //   4. Sort the chromosome's StitchedERs by start so write_to_gtf reads
    //      them in genomic order even when chains came from different passes.
    //
    // expressed_regions is an unordered_map, so its iteration order is
    // implementation-defined. Sort the keys here to make cross-chromosome
    // output order deterministic across runs and platforms.
    std::vector<std::string> chroms_sorted;
    chroms_sorted.reserve(expressed_regions.size());
    for (const auto& chrom_ers : expressed_regions) chroms_sorted.push_back(chrom_ers.first);
    std::sort(chroms_sorted.begin(), chroms_sorted.end());
    for (const std::string& chrom : chroms_sorted)
    {
        const auto& ers = expressed_regions.at(chrom);
        if (ers.empty()) continue;

        std::vector<uint32_t> plus_sjs;
        std::vector<uint32_t> minus_sjs;
        const auto sjs_it = mm_chrom_sj.find(chrom);
        if (sjs_it != mm_chrom_sj.end())
        {
            plus_sjs.reserve(sjs_it->second.size());
            minus_sjs.reserve(sjs_it->second.size());
            for (uint32_t sj_id : sjs_it->second)
            {
                if (sj_id == 0 || sj_id - 1 >= rr_all_sj.size()) continue;
                (rr_all_sj[sj_id - 1].strand ? plus_sjs : minus_sjs).emplace_back(sj_id);
            }
        }

        std::unordered_set<int> consumed;
        std::vector<StitchedER> chrom_stitched;
        if (!plus_sjs.empty())  stitch_one_strand(chrom, '+', plus_sjs,  ers, rr_all_sj, consumed, chrom_stitched);
        if (!minus_sjs.empty()) stitch_one_strand(chrom, '-', minus_sjs, ers, rr_all_sj, consumed, chrom_stitched);

        for (int i = 0; i < static_cast<int>(ers.size()); ++i)
        {
            if (consumed.count(i)) continue;
            chrom_stitched.emplace_back(StitchedER(ers[i], i));
        }

        std::sort(chrom_stitched.begin(), chrom_stitched.end(),
                  [](const StitchedER& a, const StitchedER& b) { return a.start < b.start; });
        for (auto& ser : chrom_stitched) stitched_ERs.emplace_back(std::move(ser));
    }
}


void Integrator::write_to_gtf(const std::string& output_path)
{
    std::ofstream out(output_path);
    if (!out.is_open()) {
        std::cerr << "[ERROR] could not open output file " << output_path << std::endl;
        return;
    }

    auto now = std::chrono::system_clock::now(); // get today's date
    std::chrono::year_month_day ymd{std::chrono::floor<std::chrono::days>(now)}; // formatted as YYYY-MM-DD

    // convert to string to avoid errors
    std::string date =
        std::to_string(int(ymd.year())) + "-" +
        std::to_string(unsigned(ymd.month())) + "-" +
        std::to_string(unsigned(ymd.day()));

    // write headers
    out << "#description: expressed region annotation of the genome based on BedGraph and MM / RR splice junction information." << std::endl;
    out << "#provider: FASTDER" << std::endl;
    out << "#contact: martina.lavanya@gmail.com" << std::endl;
    out << "#format: gtf" << std::endl;
    out << "#date: " << date << std::endl;

    for (unsigned int i = 0; i < this->stitched_ERs.size(); ++i)
    {
        // each stitched_er is both a gene and a transcript
        GTFRow gtf_row = GTFRow(stitched_ERs[i], "gene", i + 1);
        out << gtf_row << std::endl;
        gtf_row.change_feature("transcript", i + 1, 0);
        out << gtf_row << std::endl;
        int exon_nr = 1;
        // add the ERs within the stitched_er
        for (unsigned int k = 0; k < stitched_ERs[i].er_ids.size(); ++k)
        {
            if (stitched_ERs[i].er_ids.at(k) != -1){
                gtf_row.change_feature("exon", i + 1, exon_nr);
                // need to include SJ length as well
                gtf_row.end = gtf_row.start + stitched_ERs.at(i).all_coverages.at(k).first; // start + length = end
                gtf_row.score = stitched_ERs.at(i).all_coverages.at(k).second; // use the per-exon average coverage here instead of the overall coverage
                out << gtf_row << std::endl;
                gtf_row.start = gtf_row.end;
                ++exon_nr;
            }
            else
            {
                gtf_row.start += stitched_ERs.at(i).all_coverages.at(k).first; // add length of the SJ
            }
        }
    }
    out.close();
}
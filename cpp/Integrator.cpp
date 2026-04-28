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

// Build chains of ERs connected by SJs of one strand, and append the
// chains of length 2 or more to stitched_ERs (tagged with that strand).
// Single-ER candidates produced along the way are not emitted here; the
// caller is responsible for emitting unstitched ERs once with strand '.'.
//
// Algorithm: walk the chromosome's ERs left to right. Keep a candidate chain
// in progress. For each new ER, try to extend the candidate via the next
// available SJ on this strand. If the SJ matches (coverage and position), the
// chain grows; otherwise the candidate is closed (emitted only if it has 2+
// ERs) and a fresh candidate is seeded with the current ER. Indices of any
// ERs that wound up inside an emitted chain are added to consumed_indices.
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

    StitchedER candidate(ers[0], 0);
    candidate.strand = strand;
    bool have_candidate = true;

    auto close_candidate = [&]()
    {
        // a chain has been extended at least once when er_ids has more than the seed entry
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

    int max_chain_len = 0; // number of ERs in the longest chain emitted

    for (int i = 1; i < static_cast<int>(ers.size()); ++i)
    {
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
            // no SJs of this strand left: close the current candidate and start fresh
            // with this ER as its own (unstitched) seed
            close_candidate();
            candidate = StitchedER(expressed_region, i);
            candidate.strand = strand;
            have_candidate = true;
            continue;
        }

        // skip past SJs whose end coordinate is too far behind the chain
        while (sj_it != strand_sjs.end()
            && (candidate.end > rr_all_sj[*sj_it - 1].start
                && !within_threshold(candidate.end, rr_all_sj[*sj_it - 1].start))
            && rr_all_sj[*sj_it - 1].chrom == chrom)
        {
            ++sj_it;
        }
        // if we ran out, fall back to the last entry so dereferencing is safe
        const auto sj_to_check = (sj_it == strand_sjs.end()) ? std::prev(strand_sjs.end()) : sj_it;

        if (is_similar(candidate, expressed_region, rr_all_sj[*sj_to_check - 1]))
        {
            uint64_t sj_length = expressed_region.start - ers[candidate.er_ids.back()].end;
            candidate.append(-1, sj_length, 0.0); // append the intron placeholder
            candidate.append(i, expressed_region.length, expressed_region.coverage);
            // count actual ERs in chain (non-negative er_ids entries)
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
    // For each chromosome:
    //   1. Bucket the chromosome's SJs by strand. SJRow.strand is true for '+', false for '-'.
    //   2. Run stitch_one_strand twice (once per non-empty bucket). Each pass
    //      emits only chains of 2 or more ERs, tagged with its strand, and
    //      records the indices of consumed ERs.
    //   3. Walk the chromosome's ERs and emit each ER not consumed by any
    //      chain as a single-ER StitchedER tagged with strand '.'. This makes
    //      sure every ER appears in stitched_ERs exactly once: either inside
    //      a strand-labelled chain (when SJs supported it) or as an unstranded
    //      standalone ER (when no SJ on either strand stitched it).
    for (const auto& sjs : mm_chrom_sj)
    {
        const std::string& chrom = sjs.first;
        if (expressed_regions.find(chrom) == expressed_regions.end()) continue;
        const auto& ers = expressed_regions.at(chrom);
        if (ers.empty()) continue;

        std::vector<uint32_t> plus_sjs;
        std::vector<uint32_t> minus_sjs;
        plus_sjs.reserve(sjs.second.size());
        minus_sjs.reserve(sjs.second.size());
        for (uint32_t sj_id : sjs.second)
        {
            if (sj_id == 0 || sj_id - 1 >= rr_all_sj.size()) continue;
            (rr_all_sj[sj_id - 1].strand ? plus_sjs : minus_sjs).emplace_back(sj_id);
        }

        std::unordered_set<int> consumed;
        // Collect this chromosome's StitchedERs in a local vector so we can
        // sort them by start position before appending. The original
        // (unstranded) algorithm always emitted ERs in genomic order, and
        // write_to_gtf relies on it; running two strand passes plus a
        // standalone pass would otherwise interleave them out of order.
        std::vector<StitchedER> chrom_stitched;
        if (!plus_sjs.empty())  stitch_one_strand(chrom, '+', plus_sjs,  ers, rr_all_sj, consumed, chrom_stitched);
        if (!minus_sjs.empty()) stitch_one_strand(chrom, '-', minus_sjs, ers, rr_all_sj, consumed, chrom_stitched);

        // any ER not pulled into a chain is emitted once, unstranded
        for (int i = 0; i < static_cast<int>(ers.size()); ++i)
        {
            if (consumed.count(i)) continue;
            StitchedER fresh(ers[i], i);
            // strand stays '.' (default)
            chrom_stitched.emplace_back(fresh);
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
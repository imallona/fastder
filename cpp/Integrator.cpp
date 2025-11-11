//
// Created by martinalavanya on 20.10.25.
//

#include <cassert>
#include <Integrator.h>

// constructor
Integrator::Integrator(double coverage_tolerance_, int position_tolerance_)
{
    stitched_ERs = {}; // empty vector with elements of type StitchedER
    coverage_tolerance = coverage_tolerance_;
    position_tolerance = position_tolerance_;
}



// function that calculates relative match with a tolerance of +/- n%
bool Integrator::within_threshold(double val1, double val2){
    double tolerance_bottom = val1 * (1 - coverage_tolerance);
    double tolerance_top = val1 * (1 + coverage_tolerance);
    return val2 >= tolerance_bottom && val2 <= tolerance_top;
}


bool Integrator::within_threshold(uint64_t pos_1, uint64_t pos_2){
    //example: pos_1 = 15 (end position of exon), pos_2 = 18 (start position of SJ)
    // pos_1 must be within 18 +- 5
    // 15 >= 18 - 5 and 15 <= 18 + 5
    return pos_1 >= pos_2 - position_tolerance && pos_1 <= pos_2 + position_tolerance;
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

void Integrator::stitch_up(std::unordered_map<std::string, std::vector<BedGraphRow>>& expressed_regions, const std::map<std::string, std::vector<uint64_t>>& mm_chrom_sj, const std::vector<SJRow>& rr_all_sj)
{
    // print all splice junctions in chromosome 1
    std::cout << "--TEST SJ--" << std::endl;
    for (auto& cov : mm_chrom_sj.at("chr1"))
    {
        std::cout << cov << ", " << rr_all_sj.at(cov - 1) << std::endl;
    }
    std::cout << "--TEST SJ--" << std::endl;

    // print all ERs in chromosome 1
    std::cout << "--TEST ERs--" << std::endl;
    for (auto& er : expressed_regions.at("chr1"))
    {
        er.print();
    }
    std::cout << "--TEST ERs--" << std::endl;
    // iterate over chromosomes and sj_ids -> sjs.first = chrom, sjs.second = vector<sj_id>
    for (auto& sjs : mm_chrom_sj)
    {
        std::string chrom = sjs.first;
        StitchedER er1 = StitchedER(expressed_regions.at(chrom).at(0), 0); // define the first StitchedER, currently consisting of 1 ER
        stitched_ERs.push_back(er1);

        expressed_regions.at(chrom).at(0).print();
        auto current_sj_id = sjs.second.begin(); // iterator over the vector of sj_id

        int max_stitched_ers = 0;
        int nof_stitched_ers = 0;
        // iterate over expressed regions starting with region 2 (since region 1 was already appended)
        for (int i = 1; i < expressed_regions.at(chrom).size(); ++i)
        {
            //only compare if we aren't at the last SJ yet
            if (current_sj_id != sjs.second.end()){

                const auto& expressed_region = expressed_regions[chrom][i];
                StitchedER& current_stitched_er = stitched_ERs.back(); // this is one expressed region right now

                // skip to SJ with coordinates that line up with the most recent ER
                std::cout << "upstream ER: " << expressed_regions[chrom][current_stitched_er.er_ids.back()].chrom << ", (pos) " << expressed_regions[chrom][current_stitched_er.er_ids.back()].start << "\t" << expressed_regions[chrom][current_stitched_er.er_ids.back()].end << ", (len) " << expressed_regions[chrom][current_stitched_er.er_ids.back()].end -  expressed_regions[chrom][current_stitched_er.er_ids.back()].start <<std::endl;
                std::cout << "downstream ER: " << expressed_region.chrom << ", (pos) " << expressed_region.start << "\t" << expressed_region.end << ", (len) " <<  expressed_region.end  - expressed_region.start << std::endl;
                std::cout << "current SJ = " << rr_all_sj[*current_sj_id - 1].start << " <--> " <<  rr_all_sj[*current_sj_id - 1].end << std::endl;
                while (current_sj_id != sjs.second.end()
                    && (current_stitched_er.end > rr_all_sj[*current_sj_id - 1].start && !within_threshold(current_stitched_er.end, rr_all_sj[*current_sj_id - 1].start))
                    && rr_all_sj[*current_sj_id - 1].chrom == chrom) //& !within_threshold(most_recent_er.end, rr_all_sj[*current_sj_id - 1].start)
                {
                    std::cout << "current SJ = " << rr_all_sj[*current_sj_id - 1].start << " <--> " <<  rr_all_sj[*current_sj_id - 1].end << std::endl;
                    ++current_sj_id;

                }
                // make sure to never dereference the end() pointer
                if (current_sj_id == sjs.second.end())
                {
                    --current_sj_id;
                }

                // get rr_all_sj, which is a vector of SJRows
                if (is_similar(current_stitched_er, expressed_region, rr_all_sj[*current_sj_id - 1]))
                {

                    std::cout << "[JUNCTION] " << expressed_regions[chrom][current_stitched_er.er_ids.back()].end << " <--> " << rr_all_sj[*current_sj_id - 1].start << ", " << rr_all_sj[*current_sj_id - 1].end<< " <--> " <<  expressed_region.start << std::endl;

                    //expressed_region.print();
                    // the chromosome that
                    assert(rr_all_sj[*current_sj_id - 1].chrom == expressed_region.chrom && expressed_region.chrom == expressed_regions[chrom][current_stitched_er.er_ids.back()].chrom);
                    //expressed_region.print();
                    uint64_t sj_length = expressed_region.start - expressed_regions[chrom][current_stitched_er.er_ids.back()].end; // always use ER coordinates since a small mismatch of SJ and ER coordinates is tolerated
                    // append the spliced region and the exon
                    current_stitched_er.append(-1, sj_length, 0.0);
                    current_stitched_er.append(i, expressed_region.length, expressed_region.coverage);
                    ++nof_stitched_ers;

                    assert(current_stitched_er.end == expressed_region.end);

                    // move to next SJ
                    ++current_sj_id;

                    // find maximum number of ERs that were stitched together
                    if (max_stitched_ers < nof_stitched_ers)
                    {
                        max_stitched_ers = nof_stitched_ers;
                    }
                }

                // current ER doesn't belong to any existing ERs --> start a new ER
                else
                {
                    nof_stitched_ers = 1; // reset counter
                    stitched_ERs.push_back(StitchedER(expressed_region, i));


                }
                std::cout << "\n";


            }
        }
        std::cout << "max_stitched_ers = " << max_stitched_ers<< std::endl;
    }
}


void Integrator::write_to_gtf(const std::string& output_path)
{
    std::ofstream out(output_path);
    if (!out.is_open()) {
        std::cerr << "Error: could not open output file " << output_path << std::endl;
        return;
    }
    // get today's date
    auto now = std::chrono::system_clock::now();
    std::chrono::year_month_day ymd{std::chrono::floor<std::chrono::days>(now)};
    //format as YYYY-MM-DD
    std::string date = std::format("{:%Y-%m-%d}", ymd);

    // write headers
    out << "##description: expressed region annotation of genome based on bigwig and MM / RR splice junction information." << std::endl;
    out << "##provider: FASTDER" << std::endl;
    out << "##contact: marlehmann@ethz.ch" << std::endl;
    out << "##format: gtf" << std::endl;
    out << "##date: " << date << std::endl;

    for (unsigned int i = 0; i < this->stitched_ERs.size(); ++i)
    {
        // each stitched_er is both a gene and a transcript
        GTFRow gtf_row = GTFRow(stitched_ERs[i], "gene", i + 1);
        out << gtf_row << std::endl;

        gtf_row.change_feature("transcript", i + 1, 0);
        out << gtf_row << std::endl;
        std::cout << "gtf coords: " << gtf_row.start << " " << gtf_row.end << std::endl;
        std::cout << "stitched er coords: " << stitched_ERs[i].start << " " << stitched_ERs[i].end << std::endl;
        for (auto x : stitched_ERs[i].all_coverages)
        {
            std::cout << x.first << ", " << x.second << std::endl;
        }

        for (auto x : stitched_ERs[i].er_ids)
        {
            std::cout << x << std::endl;
        }
        // add the ERs within the stitched_er
        for (unsigned int k = 0; k < stitched_ERs[i].er_ids.size(); ++k)
        {
            if (stitched_ERs[i].er_ids.at(k) != -1){
                gtf_row.change_feature("exon", i + 1, k + 1);
                std::cout  << "start= " << gtf_row.start << ", length = " << stitched_ERs.at(i).all_coverages.at(k).first << std::endl;
                // need to use the SJ length as well
                gtf_row.end = gtf_row.start + stitched_ERs.at(i).all_coverages.at(k).first; // start + length = end
                gtf_row.score = stitched_ERs.at(i).all_coverages.at(k).second; // use the per-exon average coverage here instead of the overall coverage

                out << gtf_row << std::endl;
                gtf_row.start = gtf_row.end;
            }
            else
            {
                gtf_row.start += stitched_ERs.at(i).all_coverages.at(k).first; // add length of the SJ
                std::cout << "added sj length: " << stitched_ERs.at(i).all_coverages.at(k).first << std::endl;
            }
        }

        //gtf_row.end = stitched_ERs[i].end; // missing the last exon in the chain

    }
    out.close();
}
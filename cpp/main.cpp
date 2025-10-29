//
// Created by martinalavanya on 08/10/2025.
//

#include <iostream>

#include "BedGraphRow.h"
#include "Parser.h"
#include "Averager.h"
#include "GTFRow.h"
#include "Integrator.h"

int main() {

    // parse files
    std::string directory = "../data/all_chr";
    std::cout << "Enter directory name: ";
    //std::cin >> directory;


    // parse files
    Parser parser(directory);
    parser.search_directory();

    // get mean coverage vector
    Averager averager;
    averager.compute_mean_coverage(parser.all_per_base_coverages, parser.chromosome_sequence);

    // get expressed regions
    averager.find_ERs(0.25, 5, parser.chromosome_sequence);
    //std::cout << averager.expressed_regions.size() << std::endl;
    std::cout << " first 20 out of " << averager.expressed_regions.size() <<" expressed regions" << std::endl;
    std::cout << "chrom" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" << "\t" << "length" << std::endl;

    for (int i = 0; i  < 20; ++i)
    {
        averager.expressed_regions[parser.chromosome_sequence[0]][i].print();
    }
    std::cout << "*****************************" << std::endl;
    std::cout << " nr of sj = " << parser.mm_sj_counts.size() << std::endl;
    auto it = parser.mm_sj_counts.begin();
    for (int i = 0; i  < 20; ++i)
    {
        std::cout << parser.rr_all_sj[it->first] << std::endl;
        ++it;
    }
    // use splice junctions to stitch together expressed regions
    Integrator integrator = Integrator();
    integrator.stitch_up(averager.expressed_regions, parser.mm_sj_counts, parser.rr_all_sj, parser.chromosome_sequence);
    //
    // std::cout << "stitched ER index" << "\t" << "(" <<  "length" <<"," << "average coverage" << ")" << std::endl;
    // std::cout << "stitched_er.across_er_coverage" << "\t" << "stitched_er.start" << "\t" << "stitched_er.end" << "\t" << "stitched_er.total_length" << std::endl;
    //
    // std::cout << averager.expressed_regions.size() <<" expressed regions" << std::endl;
    // std::cout << "number of stitched regions = "  << integrator.stitched_ERs.size() << std::endl;
    // std::cout << "nr of sj in provided data + permitted chromosomes = " << parser.mm_sj_counts.size() << std::endl;
    // int chr1_count = 0;
    // int chr9_count = 0;
    // int chr19_count = 0;
    // for (auto& sj : parser.mm_sj_counts)
    // {
    //     if (parser.rr_all_sj[sj.first].chrom == "chr1")
    //     {
    //         ++chr1_count;
    //     }
    //     else if (parser.rr_all_sj[sj.first].chrom == "chr9")
    //     {
    //         ++chr9_count;
    //     }
    //     else if (parser.rr_all_sj[sj.first].chrom == "chr19")
    //     {
    //         ++chr19_count;
    //     }
    //     else
    //     {
    //         std::cout << parser.rr_all_sj[sj.first].chrom << std::endl;
    //     }
    // }
    // std::cout << "splice junctions in chr 1: " << chr1_count << std::endl;
    // std::cout << "splice junctions in chr 9: " << chr9_count << std::endl;
    // std::cout << "splice junctions in chr 19: " << chr19_count << std::endl;
    // chr1_count = 0;
    // chr9_count = 0;
    // chr19_count = 0;
    // for (auto& er : averager.expressed_regions)
    // {
    //     if (er.chrom == "chr1")
    //     {
    //         ++chr1_count;
    //     }
    //     else if (er.chrom == "chr9")
    //     {
    //         ++chr9_count;
    //     }
    //     else if (er.chrom == "chr19")
    //     {
    //         ++chr19_count;
    //     }
    // }
    // std::cout << "expressed regions in chr 1: " << chr1_count << std::endl;
    // std::cout << "expressed regions in chr 9: " << chr9_count << std::endl;
    // std::cout << "expressed regions in chr 19: " << chr19_count << std::endl;
    //
    // chr1_count = 0;
    // chr9_count = 0;
    // chr19_count = 0;
    // for (auto& stitched_er : integrator.stitched_ERs)
    // {
    //     if (stitched_er.chrom == "chr1")
    //     {
    //         ++chr1_count;
    //     }
    //     else if (stitched_er.chrom == "chr9")
    //     {
    //         ++chr9_count;
    //     }
    //     else if (stitched_er.chrom == "chr19")
    //     {
    //         ++chr19_count;
    //     }
    // }
    // std::cout << "stitched regions in chr 1: " << chr1_count << std::endl;
    // std::cout << "stitched regions in chr 9: " << chr9_count << std::endl;
    // std::cout << "stitched regions in chr 19: " << chr19_count << std::endl;
    //
    // // convert to GTF format
    // std::string output_path = "../data/output.gtf";
    // integrator.write_to_gtf(output_path);


    return 0;
}

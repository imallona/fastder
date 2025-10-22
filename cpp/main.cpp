//
// Created by martinalavanya on 08/10/2025.
//

#include <iostream>

#include "BedGraphRow.h"
#include "Parser.h"
#include "Averager.h"
#include "Integrator.h"

int main() {

    // parse files
    std::string directory = "../data/chr_19";
    std::cout << "Enter directory name: ";
    //std::cin >> directory;


    // parse files
    Parser parser(directory);
    parser.search_directory();


    Averager averager;
    // get per-base coverage
    averager.get_all_per_base_coverage(parser.all_bedgraphs);

    // get mean coverage vector
    averager.compute_mean_coverage();

    // get expressed regions
    averager.find_ERs(0.25, 5);
    //std::cout << averager.expressed_regions.size() << std::endl;
    std::cout << "chrom" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" << "\t" << "total_reads" <<  "\t" << "length" << std::endl;
    std::cout << " first 20 out of " << averager.expressed_regions.size() <<" expressed regions" << std::endl;
    for (int i = 0; i  < 20; ++i)
    {
        averager.expressed_regions[i].print();
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
    integrator.stitch_up(averager.expressed_regions, parser.mm_sj_counts, parser.rr_all_sj);

    std::cout << "stitched ER index" << "\t" << "(" <<  "length" <<"," << "average coverage" << ")" << std::endl;
    std::cout << "stitched_er.across_er_coverage" << "\t" << "stitched_er.start" << "\t" << "stitched_er.end" << "\t" << "stitched_er.total_length" << std::endl;


    for (auto & stitched_er : integrator.stitched_ERs)
    {
        std::cout << stitched_er << std::endl;
    }

    return 0;
}

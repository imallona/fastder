//
// Created by martinalavanya on 08/10/2025.
//

#include <cassert>
#include <iostream>

#include "BedGraphRow.h"
#include "Parser.h"
#include "Averager.h"
#include "GTFRow.h"
#include "Integrator.h"

int main(int argc, char* argv[]) {
    // parse command-line arguments
    std::vector<std::string> chromosomes;
    // default values (if not provided by user)
    int position_tolerance = 5;
    double coverage_tolerance = 0.1;
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "--chr")
        {
            ++i;
            while (i < argc && std::string(argv[i]).rfind("--", 0) != std::string::npos )
            {
                chromosomes.push_back(argv[i]);
                ++i;
            }
            // decrement i again
            --i;

        }
        else if (arg == "--position-tol")
        {
            position_tolerance = atoi(argv[i]);
        }

        else if (arg == "--coverage-tol")
        {
            coverage_tolerance = std::stod(argv[i]);
        }
        else
        {
            std::cout << "unknown argument '" << argv[0] << "'" << std::endl;
        }
    }
    // parse files
    std::string directory = "../data/test_exon_skipping";
    std::cout << "Enter directory name: ";
    //std::cin >> directory;


    // parse files
    Parser parser(directory, chromosomes);
    parser.search_directory();

    // get mean coverage vector
    Averager averager;
    averager.compute_mean_coverage(parser.all_per_base_coverages);
    for (const auto& cov : parser.all_per_base_coverages[0])
    {
        std::cout << cov.first << ": ";
        for (const auto& val : cov.second)
            std::cout << val << " ";
        std::cout << std::endl;
    }
    // get expressed regions
    averager.find_ERs(0.25, 5);
    //std::cout << averager.expressed_regions.size() << std::endl;
    // std::cout << " first 20 out of " << averager.expressed_regions.size() <<" expressed regions" << std::endl;
    // std::cout << "chrom" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" << "\t" << "length" << std::endl;

    for (int i = 0; i  < averager.expressed_regions.size(); ++i)
    {
        averager.expressed_regions.at(parser.chromosomes.at(0)).at(i).print();
    }
    std::cout << "*****************************" << std::endl;
    std::cout << " nr of chromosomes: " << parser.mm_chrom_sj.size() << std::endl;
    for (auto& chrom : parser.mm_chrom_sj)
    {
        std::cout << " chrom " << chrom.first << " : " << chrom.second.size() << std::endl;
    }
    // auto it = parser.mm_chrom_sj.begin();
    // for (int i = 0; i  < 20; ++i)
    // {
    //     assert(it->second.size() >= 20);
    //     std::cout << parser.rr_all_sj[it->second.at(i)] << std::endl;
    //
    // }
    // use splice junctions to stitch together expressed regions
    Integrator integrator = Integrator(coverage_tolerance, position_tolerance);
    integrator.stitch_up(averager.expressed_regions, parser.mm_chrom_sj, parser.rr_all_sj);

    // std::cout << "stitched ER index" << "\t" << "(" <<  "length" <<"," << "average coverage" << ")" << std::endl;
    // std::cout << "stitched_er.across_er_coverage" << "\t" << "stitched_er.start" << "\t" << "stitched_er.end" << "\t" << "stitched_er.total_length" << std::endl;

    std::cout << "number of stitched regions = "  << integrator.stitched_ERs.size() << std::endl;

    // SUMMARY OF SPLICE JUNCTIONS
    for (auto& c : parser.mm_chrom_sj)
    {
        std::cout << "splice junctions in chr " << c.first << " = " << c.second.size() << std::endl;
    }

    // SUMMARY OF EXPRESSED REGIONS
    for (auto& c : averager.expressed_regions)
    {
        std::cout << "expressed regions in chr " << c.first << " = " << c.second.size() << std::endl;
    }

    // SUMMARY OF STITCHED EXPRESSED REGIONS
    std::unordered_map<std::string, int> ser_counts;
    for (auto& stitched_er : integrator.stitched_ERs)
    {
        ser_counts[stitched_er.chrom]++;
    }
    for (auto& c : ser_counts)
    {
        std::cout << "stitched ERs in chr " << c.first << " = " << c.second << std::endl;
    }


    // convert to GTF format
    std::string output_path = "../data/output.gtf";
    integrator.write_to_gtf(output_path);


    return 0;
}

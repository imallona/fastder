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
#include <chrono>
int main(int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> chromosomes; //= {"chr19"}; // "chr1", "chr9",
    // default values (if not provided by user)
    int position_tolerance = 20; // [3,5,7,10]
    int min_length = 10; //[5, 10]
    double coverage_tolerance = 0.8; //[0.
    double min_coverage = 0.05;
    std::string directory  = "../data";
    bool directory_provided = true;

    std::cout
    << "\n "
    << "\t \t \t \t \t\tWELCOME TO \n"
        <<
    "        _____ _    ____ _____ ____  _____ ____  \n"
    "       |  ___/ \\  / ___|_   _|  _ \\| ____|  _ \\ \n"
    "       | |_ / _ \\ \\___ \\ | | | | | |  _| | |_) |\n"
    "       |  _/ ___ \\ ___) || | | |_| | |___|  _ < \n"
    "       |_|/_/   \\_\\____/ |_| |____/|_____|_| \\_\\"
    << "\n \n "
    << std::endl;
    std::cout << "Usage: fastder [options]\n\n"
            << "Options:\n"
            << "  --dir <path> ...             [REQUIRED] Relative path from build directory to file directory. \n"
            << "                               Example: --dir ../../data/test_exon_skipping \n\n"
            << "  --chr <chr1> <chr2> ...      List of chromosomes to process. Default = ALL\n"
            << "                               Example: --chr chr1 chr2 chr3\n\n"
            << "  --min-coverage <float> Coverage threshold to qualify as an expressed region (ER), in [CPM]. Normalization is done in-place by library size. Default = 0.25 CPM.\n"
            << "                               Example: --min-coverage 0.25\n\n"
            << "  --position-tolerance <int>   Maximum permitted position deviation of splice junction and ER coordinates, in [bp]. Default = 5 bp\n"
            << "                               Example: --position-tolerance 5\n\n"
            << "  --coverage-tolerance <float> Permitted coverage deviation between stitched ERs, as a proportion (e.g. 0.1 = 10 %). Default = 0.1\n"
            << "                               Example: --coverage-tolerance 0.1\n\n"
            << "  --help                       Show this help message.\n\n"
            << "Example:\n"
            << "  ./fastder  --dir ../../data/test_exon_skipping --chr chr1 chr2 --position-tolerance 5 "
             "--min-coverage 0.25 --coverage-tolerance 0.1\n"
            << std::endl;

    // parse command-line arguments
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "--chr") // --chr chr1 chr2
        {
            ++i;
            // the next argument is indicated with the --
            while (i < argc && std::string(argv[i]).rfind("--", 0) == std::string::npos)
            {
                chromosomes.push_back(argv[i]);
                ++i;
            }
            // decrement i again
            --i;

        }
        else if (arg == "--position-tolerance")
        {
            position_tolerance = atoi(argv[++i]);
        }
        else if (arg == "--min-length")
        {
            min_length = atoi(argv[++i]);
        }
        else if (arg == "--min-coverage")
        {
            min_coverage = std::stod(argv[++i]);
        }

        else if (arg == "--coverage-tolerance")
        {
            coverage_tolerance = std::stod(argv[++i]);
        }

        else if (arg == "--dir")
        {
            directory = argv[++i];
            directory_provided = true;
        }
        else if (arg == "--help")
        {
            std::cout << "Usage: fastder [options]\n\n"
                        << "Options:\n"
                        << "  --dir <path> ...             [REQUIRED] Relative path from build directory to file directory. \n"
                        << "                               Example: --dir ../../data/test_exon_skipping \n\n"
                        << "  --chr <chr1> <chr2> ...      List of chromosomes to process. Default = ALL\n"
                        << "                               Example: --chr chr1 chr2 chr3\n\n"
                        << "  --min-coverage <float>       Coverage threshold to qualify as an expressed region (ER), in [CPM]. Normalization is done in-place by library size. \n"
                           "                               Default = 0.25 CPM.\n"
                        << "                               Example: --min-coverage 0.25\n\n"
                        << "  --position-tolerance <int>   Maximum permitted position deviation of splice junction and ER coordinates, in [bp]. Default = 5 bp\n"
                        << "                               Example: --position-tolerance 5\n\n"
                        << "  --coverage-tolerance <float> Permitted coverage deviation within a stitched ER, as a proportion (e.g. 0.1 = 10 %). Default = 0.1\n"
                        << "                               Example: --coverage-tolerance 0.1\n\n"
                        << "  --help                       Show this help message.\n\n"
                        << "Example:\n"
                        << "  ./fastder --chr chr1 chr2 --position-tolerance 5 "
                         "--min-coverage 0.25 --coverage-tolerance 0.1\n"
                        << std::endl;

        }
        else
        {
            std::cout << "[ERROR] Unknown argument '" << argv[i] << "'" << std::endl;
        }
    }

    if (!directory_provided)
    {
        std::cout << "[ERROR] No working directory provided! Quitting...";
        return 1;
    }
    // for (auto chr : chromosomes)
    // {
    //     std::cout << chr << std::endl;
    // }

    // parse files
    Parser parser(directory, chromosomes);
    parser.search_directory();

    // get mean coverage vector
    Averager averager;
    averager.compute_mean_coverage(parser.all_per_base_coverages);
    // for (const auto& cov : parser.all_per_base_coverages[0])
    // {
    //     std::cout << cov.first << ": ";
    //     for (const auto& val : cov.second)
    //         std::cout << val << " ";
    //     std::cout << std::endl;
    // }
    // get expressed regions
    averager.find_ERs(min_coverage, min_length);
    //std::cout << averager.expressed_regions.size() << std::endl;
    // std::cout << " first 20 out of " << averager.expressed_regions.size() <<" expressed regions" << std::endl;
    // std::cout << "chrom" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" << "\t" << "length" << std::endl;

    // for (int i = 0; i  < averager.expressed_regions.size(); ++i)
    // {
    //     averager.expressed_regions.at(parser.chromosomes.at(0)).at(i).print();
    // }
    // std::cout << "*****************************" << std::endl;
    // std::cout << " nr of chromosomes: " << parser.mm_chrom_sj.size() << std::endl;
    // for (auto& chrom : parser.mm_chrom_sj)
    // {
    //     std::cout << " chrom " << chrom.first << " : " << chrom.second.size() << std::endl;
    // }

    // use splice junctions to stitch together expressed regions
    Integrator integrator = Integrator(coverage_tolerance, position_tolerance);
    integrator.stitch_up(averager.expressed_regions, parser.mm_chrom_sj, parser.rr_all_sj);

    // std::cout << "stitched ER index" << "\t" << "(" <<  "length" <<"," << "average coverage" << ")" << std::endl;
    // std::cout << "stitched_er.across_er_coverage" << "\t" << "stitched_er.start" << "\t" << "stitched_er.end" << "\t" << "stitched_er.total_length" << std::endl;


    // SUMMARY OF SPLICE JUNCTIONS
    for (auto& c : parser.mm_chrom_sj)
    {
        std::cout << "[INFO] Splice Junctions in chr" << c.first << " : " << c.second.size() << std::endl;
    }

    // SUMMARY OF EXPRESSED REGIONS
    for (auto& c : averager.expressed_regions)
    {
        std::cout << "[INFO] Expressed Regions in chr" << c.first << " : " << c.second.size() << std::endl;
    }

    // SUMMARY OF STITCHED EXPRESSED REGIONS
    std::unordered_map<std::string, int> ser_counts;
    for (auto& stitched_er : integrator.stitched_ERs)
    {
        ser_counts[stitched_er.chrom]++;
    }

    std::cout << "[INFO] Total Stitched Regions: "  << integrator.stitched_ERs.size() << std::endl;

    for (auto& c : ser_counts)
    {
        std::cout << "[INFO] Stitched ERs in chr" << c.first << " : " << c.second << std::endl;
    }


    // convert to GTF format

    // file name: FASTDER_RESULT_POS_5_COV_THR_0.1_COV_
    std::string output_path = directory + "/FASTDER_RESULT_POS_TOL_" + std::to_string(position_tolerance) + "_MIN_COV_" + std::to_string(min_coverage) + "_COV_TOL_" + std::to_string(coverage_tolerance) + "_MIN_LENGTH_" + std::to_string(min_length) + ".gtf";
    integrator.write_to_gtf(output_path);
	    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";

    return 0;
}

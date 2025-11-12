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
    double coverage_threshold = 0.25;
    std::string directory  = "../../data/test_exon_skipping";

    bool directory_provided = false;

    std::cout
    << "\n "
    << "\t \t \t WELCOME TO \n"
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
            << "  --coverage-threshold <float> Coverage threshold to qualify as an expressed region (ER), in [CPM]. Normalization is done in-place by library size. Default = 0.25 CPM.\n"
            << "                               Example: --coverage-threshold 0.25\n\n"
            << "  --position-tolerance <int>   Maximum permitted position deviation of splice junction and ER coordinates, in [bp]. Default = 5 bp\n"
            << "                               Example: --position-tolerance 5\n\n"
            << "  --coverage-tolerance <float> Permitted coverage deviation between stitched ERs, as a proportion (e.g. 0.1 = 10 %). Default = 0.1\n"
            << "                               Example: --coverage-tolerance 0.1\n\n"
            << "  --help                       Show this help message.\n\n"
            << "Example:\n"
            << "  ./fastder  --dir ../../data/test_exon_skipping --chr chr1 chr2 --position-tolerance 5 "
             "--coverage-threshold 0.25 --coverage-tolerance 0.1\n"
            << std::endl;

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
        else if (arg == "--pos-tol")
        {
            position_tolerance = atoi(argv[++i]);
        }

        else if (arg == "--cov-thr")
        {
            coverage_threshold = std::stod(argv[++i]);
        }

        else if (arg == "--cov-tol")
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
                        << "  --coverage-threshold <float> Coverage threshold to qualify as an expressed region (ER), in [CPM]. Normalization is done in-place by library size. Default = 0.25 CPM.\n"
                        << "                               Example: --coverage-threshold 0.25\n\n"
                        << "  --position-tolerance <int>   Maximum permitted position deviation of splice junction and ER coordinates, in [bp]. Default = 5 bp\n"
                        << "                               Example: --position-tolerance 5\n\n"
                        << "  --coverage-tolerance <float> Permitted coverage deviation between stitched ERs, as a proportion (e.g. 0.1 = 10 %). Default = 0.1\n"
                        << "                               Example: --coverage-tolerance 0.1\n\n"
                        << "  --help                       Show this help message.\n\n"
                        << "Example:\n"
                        << "  ./fastder --chr chr1 chr2 --position-tolerance 5 "
                         "--coverage-threshold 0.25 --coverage-tolerance 0.1\n"
                        << std::endl;

        }
        else
        {
            std::cout << "[ERROR] Unknown argument '" << argv[i] << "'" << std::endl;
        }
    }

    // if (!directory_provided)
    // {
    //     std::cout << "[ERROR] No working directory provided! Quitting...";
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

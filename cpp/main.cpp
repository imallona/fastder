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
#include <thread>
int main(int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> chromosomes;
    // default values (if not provided by user)
    int position_tolerance = 5;
    int min_length = 10;
    // coverage_tolerance gates stitching on how similar two regions' coverage
    // is. A splice junction is the real evidence that two exons belong to one
    // transcript, and the coverage of exons within a transcript varies widely,
    // so this gate defaults effectively off (any fold-difference passes) and
    // is left as a knob only for deliberate experiments.
    double coverage_tolerance = 1000.0;
    double min_coverage = 0.05;
    std::string directory;
    int cores = 10;

    std::cout
    << "\n "

    "_____ _    ____ _____ ____  _____ ____  \n"
    "|  ___/ \\  / ___|_   _|  _ \\| ____|  _ \\ \n"
    "| |_ / _ \\ \\___ \\ | | | | | |  _| | |_) |\n"
    "|  _/ ___ \\ ___) || | | |_| | |___|  _ < \n"
    "|_|/_/   \\_\\____/ |_| |____/|_____|_| \\_\\"
    << "\n\nCopyright (c) 2025 Martina Lavanya Lehmann\n"
    << std::endl;
    std::cout << "Usage: fastder [options]\n\n"
                << "Options:\n"
                << "  --dir <path> ...             [REQUIRED] Path to directory with input files (relative path to build directory or absolute path). \n"
                << "                               Example: --dir ../../data/test_exon_skipping \n\n"
                << "  --chr <chr1> <chr2> ...      List of chromosomes to process. Default = ALL.\n"
                << "                               Example: --chr chr1 chr2 or --chr 1 2 \n\n"
                << "  --min-coverage <float>       Coverage threshold to qualify as an expressed region (ER), in [CPM]. \n"
                   "                               Normalization is done in-place by library size. \n"
                   "                               Default = 0.05 CPM.\n"
                << "                               Example: --min-coverage 0.25\n\n"
                << "  --min-length <int>           Minimum length, in [nt], for a region to qualify as an expressed region (ER). Default = 10 nt.\n"
                << "                               Example: --min-length 10\n\n"
                << "  --position-tolerance <int>   Maximum permitted position deviation of splice junction and ER coordinates, in [nt]. Default = 5 nt.\n"
                << "                               Example: --position-tolerance 5\n\n"
                << "  --coverage-tolerance <float> Permitted coverage deviation when stitching, as a proportion: a\n"
                << "                               neighbour passes if its coverage is within (1 +/- tolerance)\n"
                << "                               of the running average. Not bound to (0,1); 2.0 allows a 3x\n"
                << "                               spread. Default = 1000 (the gate is effectively off, stitching\n"
                << "                               is driven by splice junctions).\n"
                << "                               Example: --coverage-tolerance 2.0\n\n"
                << "  --cores <int>                Number of cores that fastder may use. Default = 10 cores.\n"
                << "                               Example: --cores 23\n\n"
                << "Example:\n"
                << "  ./fastder --dir ../data --chr chr1 chr2 --position-tolerance 5 "
                 "--min-coverage 0.05 --coverage-tolerance 0.7 --cores 23\n"
                << std::endl;

    // parse command-line arguments
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "--chr") // --chr chr1 chr2
        {
            ++i;
            // the next argument starts with --
            while (i < argc && std::string(argv[i]).rfind("--", 0) == std::string::npos)
            {
                std::string chrom = argv[i];
                if (chrom.find("chr") == std::string::npos)
                {
                    chrom = "chr" + chrom; // add "chr" prefix if user passes --chr 1 2 3
                }
                chromosomes.push_back(chrom);
                ++i;
            }
            // decrement i again
            --i;

        }
        else if (arg == "--position-tolerance")
        {
            position_tolerance = atoi(argv[++i]);
        }
        else if (arg == "--cores")
        {
            cores = atoi(argv[++i]);
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
        }
        else
        {
            std::cout << "[ERROR] Unknown argument '" << argv[i] << "'" << std::endl;
            return 1;
        }
    }

    // exit if no directory is provided
    if (directory.empty())
    {
        std::cerr << "[ERROR] No directory with input files specified! Quitting..." << std::endl;
        return 1;
    }

    // parse files
    std::cout << "[INFO] Expecting to parse MM, RR, BedGraph and Metadata CSV files from " << directory << std::endl;
    Parser parser(directory, chromosomes, cores);
    parser.search_directory();

    // print parsing duration
    auto end_parsing = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parsing = end_parsing - start;
    std::cout << "[INFO] Parsing took " << elapsed_parsing.count() << " seconds.\n \n";

    // get mean coverage as sparse intervals (no dense per-base expansion)
    Averager averager(cores);
    averager.compute_mean_coverage(parser.all_bedgraphs);

    averager.find_ERs(min_coverage, min_length);

    // use splice junctions to stitch together expressed regions
    Integrator integrator = Integrator(coverage_tolerance, position_tolerance);
    integrator.stitch_up(averager.expressed_regions, parser.mm_chrom_sj, parser.rr_all_sj);


    // SUMMARY OF SPLICE JUNCTIONS
    for (auto& c : parser.mm_chrom_sj)
    {
        std::cout << "[INFO] Splice junctions in " << c.first << " : " << c.second.size() << std::endl;
    }

    // SUMMARY OF EXPRESSED REGIONS
    for (auto& c : averager.expressed_regions)
    {
        std::cout << "[INFO] Expressed Regions in " << c.first << ": " << c.second.size() << std::endl;
    }

    // SUMMARY OF STITCHED EXPRESSED REGIONS
    std::unordered_map<std::string, int> ser_counts;
    for (auto& stitched_er : integrator.stitched_ERs)
    {
        ser_counts[stitched_er.chrom]++;
    }

    std::cout << "[INFO] Total Stitched ERs across all chromosomes: "  << integrator.stitched_ERs.size() << std::endl;

    for (auto& c : ser_counts)
    {
        std::cout << "[INFO] Stitched ERs in " << c.first << ": " << c.second << std::endl;
    }


    // convert to GTF format
    std::string prefix = (directory.back() == '/') ? "FASTDER_RESULT_POS_TOL_" : "/FASTDER_RESULT_POS_TOL_";
    std::string output_path = directory + prefix + std::to_string(position_tolerance) + "_MIN_COV_" + std::to_string(min_coverage) + "_COV_TOL_" + std::to_string(coverage_tolerance) + "_MIN_LENGTH_" + std::to_string(min_length) + ".gtf";
    integrator.write_to_gtf(output_path);

    // print duration
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "[INFO] Total duration (wallclock time): " << elapsed.count() << " seconds.\n";

    return 0;
}

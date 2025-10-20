//
// Created by martinalavanya on 08/10/2025.
//

#include <iostream>

#include "BedGraphRow.h"
#include "Parser.h"
#include "Averager.h"

int main() {

    // parse files
    std::string directory = "../data";
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
    std::cout << averager.expressed_regions.size() << std::endl;

    return 0;
}

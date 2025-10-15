//
// Created by martinalavanya on 08/10/2025.
//

#include <iostream>

#include "BedGraphRow.h"
#include "Parser.h"

int main() {

    // parse files
    std::string directory = "../data";
    std::cout << "Enter directory name: ";
    //std::cin >> directory;


    // parse files
    Parser parser(directory);
    parser.search_directory();

    for (auto p : parser.rail_id_mm_id)
    {
        std::cout << "rail id = " << p.first << " , MM id = " << p.second << std::endl;
    }

    return 0;
}

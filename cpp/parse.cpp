//
// Created by marti on 08/10/2025.
//
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <BedGraphRow.h>

#include "parse.h"

int read_file(const std::string filename,
    std::vector<double>& per_base_coverage,
    std::vector<BedGraphRow>& data)
{

    std::cout << filename << std::endl;
    //read in file from path
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file " << filename << std::endl;
    }
    std::string line;
    int library_size = 0;
    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);
        BedGraphRow row;
        iss >> row.chrom >> row.start >> row.end >> row.coverage;
        // calculate total number of reads that map to this bp interval
        row.total_reads += (row.end - row.start) * (row.coverage); //if start = 22, end = 25, coverage = 3 --> (25 - 22) * 3 = 3 * 3 = 9
        library_size += row.total_reads;
        //row.print();
        compute_per_base_coverage(row, per_base_coverage);
        data.push_back(row);
    }

    //return library size
    return library_size;

}
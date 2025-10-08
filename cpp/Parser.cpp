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
#include <cassert>

#include "Parser.h"

// constructor
Parser::Parser(std::string _path) {
    path = _path;
}

// parse relevant chromosomes of a bedgraph file
std::vector<BedGraphRow> Parser::read_bedgraph(const std::string filename)
{
    std::vector<BedGraphRow> bedgraph;
    std::cout << filename << std::endl;
    //read in file from path
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening .bedgraph file " << filename << std::endl;
    }
    std::string line;
    unsigned int library_size = 0;
    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);
        BedGraphRow row;
        iss >> row.chrom >> row.start >> row.end >> row.coverage;
        // calculate total number of reads that map to this bp interval
        // TODO think about int -> unsigned int type safety
        row.total_reads += (row.end - row.start) * (row.coverage); //if start = 22, end = 25, coverage = 3 --> (25 - 22) * 3 = 3 * 3 = 9
        library_size += row.total_reads;
        row.print();
        // compute_per_base_coverage(row, per_base_coverage);
        bedgraph.push_back(row);
    }

    return bedgraph;


}


// read rr file
void Parser::read_rr(std::string filename)
{
    std::cout << filename << std::endl;
    //read in file from path
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file " << filename << std::endl;
    }
    std::string line;

    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);

        // skip invalid lines, headers, ERCC and Y-chromosome
        if (line.empty() || line.find("chromosome") != std::string::npos || line.find("ERCC-") != std::string::npos || line.find("chrY") != std::string::npos) {
            std::cout << line << std::endl;
            continue;
        }
        SJRow row = SJRow();
        iss >> row;

        //std::cout << row << std::endl;

        sample_rr.push_back(row);
    }


}


// pass dictionary mm_by_samples by reference and fill it with keys (sample id) and int-int pairs (splice junction id, count)
void Parser::read_mm(std::string filename) {

        std::cout << filename << std::endl;
        //read in file from path
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error opening file " << filename << std::endl;
        }
        std::string line;
        bool seen_header = false;
        while (std::getline(file, line))
        {
            // read in line by line
            std::istringstream iss(line);


            //skip comments
            if (line[0] == '%') continue;

            // allows skipping the first line without a %
            if (!seen_header) {
                seen_header = true;
                continue;
            }

            // skip invalid lines
            unsigned int sj_id, sample_id, count;
            if (!(iss >> sj_id >> sample_id >> count)){
                std::cout << "malformed line in MM file: " << line << std::endl;
                continue;
            }
            //std::cout << sj_id << " " << sample_id << " " << count << std::endl;
            // add line to MM dictionary
            mm_by_samples[sample_id].push_back(std::make_pair(sj_id, count));
            std::cout << sample_id << ": " << mm_by_samples[sample_id][0].first << " " << mm_by_samples[sample_id][0].second << std::endl;

        }

    }

void Parser::search_directory() {

    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::string filename = entry.path().string();
        std::cout << filename << std::endl;
        // read RR file
        if (filename.find("ALL.RR") != std::string::npos) {
            std::cout << "M" << std::endl;
            read_rr(filename);

        }

        else if (filename.find("ALL.MM") != std::string::npos) {
            std::cout << "R";
            read_mm(filename);
        }

        else if (filename.find(".bedGraph") != std::string::npos) {
            std::cout << "B";
            std::vector<double> per_base_coverage;
            std::vector<BedGraphRow> bedgraph = read_bedgraph(filename);
            // add to matrix
            all_bedgraphs.push_back(bedgraph);
            //all_per_base_coverages.push_back(per_base_coverage);

        }

        else {
            std::cout << "UNKNOWN FILE CATEGORY " << filename  << std::endl;
        }

        //int library_size = read_file(entry.path().string(), per_base_coverage, bedgraph);

        // normalize read counts
        //normalize(per_base_coverage, library_size);



    }

}

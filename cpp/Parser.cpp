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
#include <algorithm>

#include "Parser.h"

#include <filesystem>

// constructor
Parser::Parser(std::string _path) {
    path = _path;
}


// parse relevant chromosomes of a bedgraph file
std::vector<BedGraphRow> Parser::read_bedgraph(const std::string filename, unsigned int& library_size)
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
        //row.print();
        // compute_per_base_coverage(row, per_base_coverage);
        bedgraph.push_back(row);
    }

    // add the sample and its mm index (= the rank of the rail id across the study) to rail_id_to_mm
    // [&] references all necessary variables i.e. the required context, here it's filename
    auto it = std::find_if(all_sample_ids.begin(), all_sample_ids.end(), [&](const auto& sample)
    {
        // the external id is part of the filename for all three sources GTEX, TCGA and SRA
        return filename.find(sample.second) != std::string::npos;
    });

    unsigned int mm_id = std::distance(all_sample_ids.begin(), it);
    unsigned int rail_id = it->first;

    rail_id_to_mm.push_back(std::make_pair(rail_id, mm_id));

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
            //std::cout << line << std::endl;
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
        unsigned int max_index = 0;
        unsigned int min_index = 200;
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
            if (sample_id > max_index)
            {
                max_index = sample_id;
            }
            if (sample_id < min_index)
            {
                min_index = sample_id;
            }
            if (!(iss >> sj_id >> sample_id >> count)){
                //std::cout << "malformed line in MM file: " << line << std::endl;
                continue;
            }
            //std::cout << sj_id << " " << sample_id << " " << count << std::endl;
            // add line to MM dictionary
            mm_by_samples[sample_id].push_back(std::make_pair(sj_id, count));
            //std::cout << sample_id << ": " << mm_by_samples[sample_id][0].first << " " << mm_by_samples[sample_id][0].second << std::endl;
            //if sample_id in rail_id_to_mm ->

        }
        std::cout <<  "max index" << max_index << std::endl;
        std::cout <<  "min index" << min_index << std::endl;
        std::cout << "nr of samples" << mm_by_samples.size() << std::endl;

    }

// parse bigwig URL list csv file
void Parser::read_url_csv(std::string filename)
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

        // invalid or header line
        if (line.empty() || line == "rail_id,external_id,study,BigWigURL")
        {
            continue;
        }
        //int rail_id;
        std::string rail_id_str, sample_id;
        //iss >> sample_id >> rail_id; // only read in the first two tab-separated entries, ignore the rest!
        if (std::getline(iss, rail_id_str, ',') && std::getline(iss, sample_id, ','))
        {
            //std::cout << rail_id_str << std::endl;
            int rail_id = std::stoi(rail_id_str);
            all_sample_ids.push_back(std::make_pair(rail_id, sample_id));
        }


    }

    // sorting is n log n and finding the position + inserting can be n*n, so better to push_back and then sort

    // sort all_sample_ids by rail_id to obtain the index used in the MM file
    std::sort(all_sample_ids.begin(), all_sample_ids.end(), [](const auto& a, const auto& b)
    {
        return a.first < b.first;
    });


}


// attempt to parse all files in path (not recursive!)
void Parser::search_directory() {

    bool contains_ids = false;

    // first check for the external_id to rail_id mapping CSV file
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::string filename = entry.path().string();

        if (filename.find("BigWig_list") != std::string::npos && filename.find(".csv") != std::string::npos) //TODO I checked some filenames of the URL csv files manually and they all contain the substring BigWig_list, so I hope that this is a general rule
        {
            read_url_csv(filename);
            contains_ids = true;
            break;
        }
    }
    if (!contains_ids)
    {
        std::cerr << "MISSING BigWig URL list! Cannot proceed";
        return;
    }
    // now parse all other files
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::string filename = entry.path().string();

        //std::cout << filename << std::endl;
        // read RR file
        if (filename.find("ALL.RR") != std::string::npos) {
            std::cout << "R" << std::endl;
            read_rr(filename);

        }

        else if (filename.find("ALL.MM") != std::string::npos) {
            std::cout << "M"<< std::endl;
            read_mm(filename);
        }

        else if (filename.find(".bedGraph") != std::string::npos) {
            std::cout << "B"<< std::endl;

            std::vector<double> per_base_coverage;
            unsigned int library_size = 0;
            std::vector<BedGraphRow> sample_bedgraph = read_bedgraph(filename, library_size);

            //normalize to CPM
            for (BedGraphRow row : sample_bedgraph)
            {
                row.normalize(library_size);
            }

            // add to matrix of all bedgraphs per sample
            // all_bedgraphs[sample_id].push_back(sample_bedgraph); T

        }
        // list with the mapping of external_id to rail_id
        else if (filename.find("BigWig_list") != std::string::npos && filename.find(".csv") != std::string::npos)
        {
            continue;
        }
        else {
            std::cout << "UNKNOWN FILE CATEGORY " << filename  << std::endl;
        }

        //int library_size = read_file(entry.path().string(), per_base_coverage, bedgraph);

        // normalize read counts
        //normalize(per_base_coverage, library_size);
    }
    std::cout << "DONE" << std::endl;

}
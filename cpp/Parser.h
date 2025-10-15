//
// Created by marti on 08/10/2025.
//

#ifndef FASTDER_PARSE_H
#define FASTDER_PARSE_H
#include "SJRow.h"
#include <unordered_map>


class Parser {
public:
    Parser(std::string _path);
    void search_directory();
    std::vector<BedGraphRow> read_bedgraph(std::string filename, unsigned int& library_size);
    void read_mm(std::string filename);
    void read_rr(std::string filename);
    void normalize(const unsigned int& library_size);
    void read_url_csv(std::string filename);



    std::string path;
    std::unordered_map<std::string, std::vector<BedGraphRow>> all_bedgraphs; //key = sample id, value = bedgraph of the sample
    //std::vector<std::unordered_map<std::string, std::vector<double>>> all_per_base_coverages;


    // store RR file for one sample
    std::vector<SJRow> sample_rr;

    // store Market Matrix (MM) file for one sample
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> mm_by_samples;

    std::vector<std::pair<unsigned int, unsigned int>> mm_by_sample;

    std::vector<std::pair<unsigned int, std::string>> all_sample_ids; // rail_id, sample_id

    std::vector<std::pair<unsigned int, unsigned int>> rail_id_to_mm; // rail_id, mm_id





};


#endif //FASTDER_PARSE_H
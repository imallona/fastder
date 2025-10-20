//
// Created by marti on 08/10/2025.
//

#ifndef FASTDER_PARSE_H
#define FASTDER_PARSE_H
#include "SJRow.h"
#include <unordered_map>
#include <map>


class Parser {
public:
    Parser(std::string _path);
    void search_directory();
    std::vector<BedGraphRow> read_bedgraph(const std::string& filename, uint64_t& library_size);
    void read_mm(std::string filename);
    void read_rr(std::string filename);
    void normalize(const unsigned int& library_size);
    void read_url_csv(std::string filename);
    void fill_up(std::vector<std::string> bedgraph_files);
    //void get_per_base_coverages();
    // TODO add function get_rail_id_from_filename(filename)?

    std::string path;
    std::vector<std::vector<BedGraphRow>> all_bedgraphs; //TODO maybe change to unordered map with key = sample id, value = bedgraph of the sample?

    // store RR file for all splice junctions
    std::vector<SJRow> rr_all_sj;

    // store Market Matrix (MM) file for one sample
    // std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> mm_by_samples;

    //TODO use uint64_t instead for the sj_id
    std::unordered_map<unsigned int, unsigned int> mm_sj_counts; // <sj_id, count> map for all samples in the user input

    std::vector<std::pair<unsigned int, std::string>> rail_id_to_ext_id; // <rail_id, external_id> for all samples in the dataset

    std::vector<std::pair<unsigned int, unsigned int>> rail_id_to_mm_id; // <rail_id, mm_id> mapping

};


#endif //FASTDER_PARSE_H
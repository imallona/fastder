//
// Created by marti on 08/10/2025.
//

#ifndef FASTDER_PARSE_H
#define FASTDER_PARSE_H
#include "SJRow.h"
#include <unordered_map>
#include <map>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

class Parser {
public:
    Parser(std::string _path, std::vector<std::string> chromosomes_);
    void search_directory();
    // Cache MM file because it takes so long to parse
    void save_mm_cache_(const std::filesystem::path& cache) const;
    bool load_mm_cache_(const std::filesystem::path& cache);
    void read_mm_cached_always(const std::string& filename);

    // read in individual file types
    std::vector<BedGraphRow> read_bedgraph(const std::string& filename, uint64_t& library_size);
    void read_mm(std::string filename);
    void read_rr(std::string filename);
    void read_url_csv(std::string filename);
    void fill_up(std::vector<std::string> bedgraph_files);
    [[nodiscard]] bool chr_permitted(const std::string& chr) const;

    static void compute_per_base_coverage(const BedGraphRow& row, std::unordered_map<std::string, std::vector<double>>& per_base_coverage);

    // TODO add function get_rail_id_from_filename(filename)?

    std::string path;
    std::vector<std::string> chromosomes;
    std::vector<std::vector<BedGraphRow>> all_bedgraphs; //TODO maybe change to unordered map with key = sample id, value = bedgraph of the sample?
    std::vector<std::unordered_map<std::string, std::vector<double>>> all_per_base_coverages; //NOT ordered by chromosomes
    // store RR info for each splice junctions
    std::vector<SJRow> rr_all_sj;

    // store relevant Market Matrix (MM)
    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj; // <chrom, sj_id> ordered by sj_id, map of sj occurring in samples part of the user input

    std::vector<std::pair<unsigned int, std::string>> rail_id_to_ext_id; // <rail_id, external_id> for all samples in the dataset

    std::vector<std::pair<unsigned int, unsigned int>> rail_id_to_mm_id; // <rail_id, mm_id> mapping


    const std::vector<std::string> permitted_chromosomes =  {
        "chr1",
         "chr2",
         "chr3",
         "chr4",
         "chr5",
         "chr6",
         "chr7",
         "chr8",
        "chr9",
         "chr10",
         "chr11",
         "chr12",
         "chr13",
         "chr14",
         "chr15",
         "chr16",
         "chr17",
         "chr18",
         "chr19",
        "chr20",
         "chr21",
         "chr22",
         "chrX",
    };



};


#endif //FASTDER_PARSE_H
//
// Created by marti on 08/10/2025.
//

#ifndef FASTDER_PARSE_H
#define FASTDER_PARSE_H
#include "SJRow.h"
#include <unordered_map>
#include <filesystem>
#include <fstream>
#include <unordered_set>

namespace fs = std::filesystem;

class Parser {
public:
    Parser(std::string path_, std::vector<std::string> chromosomes_, int cores_);
    void search_directory();

    // read in individual file types
    void read_all_bedgraphs(std::vector<std::string> bedgraph_files, unsigned int nof_threads);
    std::vector<BedGraphRow> read_bedgraph(const std::string& filename, uint64_t& library_size) const;
    // Parse a BigWig file into the same sparse interval representation we use
    // for BedGraph. _strand is stamped on every emitted row ('.' = unstranded,
    // '+' / '-' for stranded BigWigs). Implementation is gated on the CMake
    // option FASTDER_USE_LIBBIGWIG. Without that flag the function logs an
    // error and returns an empty vector, which keeps the build hermetic.
    std::vector<BedGraphRow> read_bigwig(const std::string& filename, uint64_t& library_size,
                                         char strand = '.') const;
    void read_mm(std::string filename);
    void read_rr(std::string filename);
    void read_url_csv(std::string filename);
    void fill_up(std::vector<std::string> bedgraph_files);

    // TODO add function get_rail_id_from_filename(filename)?
    unsigned int user_cores;
    // Minimum read support a junction needs, summed over the loaded samples,
    // to be kept for stitching. 0 keeps every junction the MM file lists,
    // which is what fastder did before the option existed.
    unsigned int min_junction_reads = 0;
    std::string path;
    std::vector<std::string> chromosomes_vec; // for fast iteration
    std::unordered_set<std::string> chromosomes_set; // for fast check if chromosome is included
    // Per-sample sparse interval coverage. Element s is the BedGraph or BigWig
    // content of sample s, sorted by (chrom, start). This is the only coverage
    // representation fastder keeps in memory; the dense per-base expansion
    // that the previous version computed has been removed since
    // compute_mean_coverage and find_ERs now operate on intervals directly.
    std::vector<std::vector<BedGraphRow>> all_bedgraphs;
    // RR rows on chromosomes the user requested. Dropped rows are not stored.
    std::vector<SJRow> rr_all_sj;

    // Maps the 1-based sj_id used by the on-disk MM file, which equals the
    // row number in the RR file, to the 1-based index into rr_all_sj. Dropped
    // rows have value 0. Size equals the total RR row count.
    std::vector<uint32_t> sj_id_remap;

    // Total junction rows seen in the RR file, independent of the chr filter.
    // Used to validate the MM header's nr_of_sj.
    uint64_t rr_total_rows = 0;

    // Sparse junction matrix per chromosome. sj_ids stored here are
    // post-remap indexes into rr_all_sj, not original RR row numbers.
    // uint32_t fits any realistic junction catalog (recount3 hg38 has 9.5M).
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;

    std::vector<std::pair<unsigned int, std::string>> rail_id_to_ext_id; // <rail_id, external_id> for all samples in the dataset
    // later sorted by rail_id to receive rank (= mm_id)

    std::unordered_set<unsigned int> mm_ids; // unordered map for fast mm_id lookup


    const std::unordered_set<std::string> permitted_chromosomes =  {
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
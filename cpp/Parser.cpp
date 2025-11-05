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
#include <cstdint> // for library size which can be too large for unsigned int

#include "Parser.h"

#include <filesystem>

// constructor
Parser::Parser(std::string _path) {
    path = _path;
}

// function to prevent reading in Y chromosome and quality check chromosomes like ERCC
bool Parser::chr_permitted(const std::string& chr) const
{
    for (auto& permitted_chromosome : this->permitted_chromosomes)
    {
        if (chr == permitted_chromosome)
        {
            return true;
        }
    }

    return false;
}
// fill vector with coverage value per bp (since different bedgraphs have different binning intervals)
void Parser::compute_per_base_coverage(const BedGraphRow& row, std::unordered_map<std::string, std::vector<double>>& per_base_coverage)
{
    // row.end is NOT inclusive
    unsigned int position = row.end - row.start; //for just one nt, row.start = 22, row.end = 23 -> position = 1
    do
    {
        // if (row.length < 10)
        //     std::cout << "row coverage = " << row.coverage << std::endl;
        per_base_coverage[row.chrom].push_back(row.coverage);
        position--;
    }
    while (position > 0);

}

// parse relevant chromosomes of a bedgraph file
std::vector<BedGraphRow> Parser::read_bedgraph(const std::string& filename, uint64_t& library_size)
{
    std::vector<BedGraphRow> bedgraph; // stores the full bedgraph of one sample, organized by rows (bins) with the same coverage
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
        BedGraphRow row = BedGraphRow();
        iss >> row.chrom >> row.start >> row.end >> row.coverage;
        if (chr_permitted(row.chrom)){
            // calculate total number of reads that map to this bp interval
            row.length = row.end - row.start;
            // end is not inclusive, since row1.end == row2.start of the next row
            row.total_reads = row.length * row.coverage; //if start = 22, end = 25, coverage = 3 --> (25 - 22) * 3 = 3 * 3 = 9
            library_size += row.total_reads;
            bedgraph.push_back(row);
        }

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

        // skip invalid lines and headers (which contain the string "chromosome" --> actually don't skip ERCC and Y-chromosome since we need all sj_ids to be correct
        if (line.empty() || line.find("chromosome") != std::string::npos) {//|| line.find("ERCC-") != std::string::npos || line.find("chrY") != std::string::npos) {
            //std::cout << line << std::endl;
            continue;
        }
        SJRow row = SJRow();
        iss >> row;

        // rr_all_sj needs to contain all sj_ids, even those of chromosomes that aren't provided in the bedgraph files --> otherwise the mapping from RR to MM file via sj_id is broken
        rr_all_sj.push_back(row);

        // store the sequence of chromosomes in the RR file --> append chromosome to vector if it's not an element of the vector yet
        // chromosome_sequence will serve as the key iteration sequence for ALL unordered maps where key = chromosome!
        if (this->chr_permitted(row.chrom) && std::ranges::find(chromosome_sequence, row.chrom) == chromosome_sequence.end()) {
            chromosome_sequence.push_back(row.chrom);
        }

    }
    std::cout << "nr of splice junctions in this study: " << rr_all_sj.size() << std::endl;
    assert(rr_all_sj.size() == 9484210);



}


// create dictionary mm_sj_counts with keys (sj_id) and values (cumulative count of this sj across samples)
// IMPORTANT: the RR file is NOT sorted by chromosomes!
void Parser::read_mm(std::string filename) {

        std::cout << filename << std::endl;
        //read in file from path
        std::ifstream file(filename);
        // max index is 2931 (= nr of samples)
        // min index is 0
        // mm_by_samples.size() = 2931

        if (!file.is_open())
        {
            std::cerr << "Error opening file " << filename << std::endl;
        }
        std::string line;
        bool seen_header = false;
        uint64_t nr_of_sj, sj_occ_in_samples, nr_of_samples;
        //auto sj_id_prev = 0;
        uint64_t count_lines = 0;
        while (std::getline(file, line))
        {


            ++count_lines;
            // read in line by line
            std::istringstream iss(line);


            //skip comments
            if (line[0] == '%') continue;

            // allows skipping the first line without a %
            if (!seen_header) {
                // header: 9484210	2931	699368828, actual #lines = 699368831
                iss >> nr_of_sj >> nr_of_samples >> sj_occ_in_samples;
                assert(nr_of_sj == rr_all_sj.size());
                seen_header = true;
                continue;
            }

            // skip invalid lines
            uint64_t sj_id;
            unsigned int mm_id, count;

            if (!(iss >> sj_id >> mm_id >> count)){
                std::cout << "malformed line in MM file: " << line << std::endl;
                std::cout << line << std::endl;
                continue;
            }

            // OLD: mm_by_samples[sample_id].push_back(std::make_pair(sj_id, count));
            // NEW: cumulative counts of a sj_id across all samples in the input

            // find the rail_id based on the mm_id --> only add sj_id if the mm_id is part of the samples
            auto it = std::find_if(rail_id_to_mm_id.begin(), rail_id_to_mm_id.end(), [&] (const auto& p)
            {
                return p.second == mm_id;
            });

            // add count if the mm was found and if chr is in bedgraph
            //std::cout << rr_all_sj[sj_id].chrom << " for splice junction " << sj_id << std::endl;

            if (it != rail_id_to_mm_id.end() && this->chr_permitted(rr_all_sj[sj_id].chrom)) // rail_id_to_mm_id has <rail_id, mm_id> mapping
            {
                // store vector of sj_ids for each chromosome
                mm_chrom_sj[rr_all_sj[sj_id].chrom].push_back(sj_id); // this creates the binding if it doesn't exist yet, initializes it to 0 and then increases it by count

            }
        }
        std::cout << "nr of lines read in MM file: " << count_lines << std::endl;
        assert(sj_occ_in_samples <= count_lines);
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
    bool first_line = true; // in case the header is parsed differently
    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);

        // invalid or header line
        if (line.empty() || line == "rail_id,external_id,study,BigWigURL" || first_line)
        {
            first_line = false;
            continue;
        }
        //int rail_id;
        std::string rail_id_str, sample_id;
        //iss >> sample_id >> rail_id; // only read in the first two tab-separated entries, ignore the rest!
        if (std::getline(iss, rail_id_str, ',') && std::getline(iss, sample_id, ','))
        {
            //convert to integer
            int rail_id = std::stoi(rail_id_str);
            rail_id_to_ext_id.push_back(std::make_pair(rail_id, sample_id));
        }


    }

    // sorting is n log n and finding the position + inserting can be n*n, so better to push_back and then sort


}

// creates a map of rail_id to mm_id in rail_id_to_mm_id
// bedgraph_files contains the file names of all samples
void Parser::fill_up(std::vector<std::string> bedgraph_files)
{
    //fill up rail_id_to_mm_id
    for (auto& bedgraph_file : bedgraph_files)
    {
        // add the sample and its mm_id (= the rank of the rail id across the study, so all files in total) to rail_id_to_mm
        // [&] references all necessary variables i.e. the required context, here it's filename
        auto it = std::find_if(rail_id_to_ext_id.begin(), rail_id_to_ext_id.end(), [&](const auto& sample)
        {
            // search for the external_id in rail_id_to_ext_id and then obtain the rail_id
            // the external id is part of the filename for all three sources GTEX, TCGA and SRA
            return bedgraph_file.find(sample.second) != std::string::npos;
        });
        if (it != rail_id_to_ext_id.end())
        {
            unsigned int mm_id = std::distance(rail_id_to_ext_id.begin(), it) + 1; // std::distance counts the steps between two iterators --> mm_id is 1 too small, so add 1
            unsigned int rail_id = it->first;

            rail_id_to_mm_id.push_back(std::make_pair(rail_id, mm_id));
        }
        else
        {
            std::cout << "ERROR: file " << bedgraph_file << "has no rail_id! " << std::endl;
        }
    }
}



// attempt to parse all files in path (not recursive!)
void Parser::search_directory() {

    bool contains_ids = false;

    // first check for the external_id to rail_id mapping CSV file
    std::vector<std::string> bedgraph_files;
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::string filename = entry.path().string();
        // create rail_id_to_ext_id
        if (filename.find("BigWig_list") != std::string::npos && filename.find(".csv") != std::string::npos) //TODO I checked some filenames of the URL csv files manually and they all contain the substring BigWig_list, so I hope that this is a general rule
        {
            std::cout << "BigWig URL list" << std::endl;
            read_url_csv(filename);
            contains_ids = true;
        }
        // read RR file
        else if (filename.find("ALL.RR") != std::string::npos) {
            std::cout << "RR file" << std::endl;
            read_rr(filename);

        }

        // collect all bedgraph files to later fill up rail_id_to_mm_id
        if (filename.find(".bedGraph") != std::string::npos)
        {
            bedgraph_files.push_back(filename);
        }
    }

    // program cannot run with missing BigWig URL list
    if (!contains_ids)
    {
        std::cerr << "MISSING BigWig URL list! Cannot proceed...";
        return;
    }

    // sort rail_id_to_ext_id by rail_id to obtain the index used in the MM file
    std::sort(rail_id_to_ext_id.begin(), rail_id_to_ext_id.end(), [](const auto& a, const auto& b)
    {
        return a.first < b.first;
    });

    std::cout << "total nr of samples in this study:  " << rail_id_to_ext_id.size() << std::endl;

    // fill up rail_id_to_mm_id mapping for all rail_ids provided by the user
    fill_up(bedgraph_files);
    std::cout << "nr of samples provided by user: " << rail_id_to_mm_id.size() << std::endl;

    // now parse all other files
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::string filename = entry.path().string();

        //std::cout << filename << std::endl;

        // don't read in cached MM files as regular MM files!
        if (entry.path().extension().string() == ".MM" && filename.find("ALL.MM") != std::string::npos && filename.find("mmcache") == std::string::npos ) {
            std::cout << "MM file"<< std::endl;
            // TODO change!
            //read_mm(filename);
            read_mm_cached_always(filename);

        }

        else if (filename.find(".bedGraph") != std::string::npos) {
            std::cout << "Bedgraph file"<< std::endl;
            //
            uint64_t library_size = 0; // ensure that the integer type is large enough
            std::vector<BedGraphRow> sample_bedgraph = read_bedgraph(filename, library_size);
            std::unordered_map<std::string, std::vector<double>> per_base_coverage;


            //normalize to CPM and expand rows to per base coverage (also normalized)
            for (BedGraphRow& row : sample_bedgraph)
            {
                row.normalize(library_size);
                compute_per_base_coverage(row, per_base_coverage);
            }

            // add to matrix of all bedgraphs per sample
            all_bedgraphs.push_back(sample_bedgraph);
            all_per_base_coverages.push_back(per_base_coverage);

        }
        else if ((filename.find("BigWig_list") != std::string::npos && filename.find(".csv") != std::string::npos) || (filename.find("ALL.RR") != std::string::npos)){
            continue;
        }
        else {
            std::cout << "UNKNOWN FILE CATEGORY " << filename  << std::endl;
        }

    }

    std::cout << "FINISHED PARSING" << std::endl;

}
// CACHING MM FILE

namespace fs = std::filesystem;

//choose caching path
static fs::path mm_cache_path(const fs::path& src) {
    //return mm_path.parent_path() / (mm_path.filename().string() + ".mmcache");
    if (src.extension() == ".mmcache") return src;
    fs::path cache = src;
    cache += ".mmcache";

    return cache;
}

//serialize parsed result (mm_chrom_sj)
void Parser::save_mm_cache_(const fs::path& cache) const {
    std::ofstream out(cache, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open cache for write: " + cache.string());

    //magic + version (helps future-proofing)
    const uint32_t magic = 0x4D4D4348; // "MMCH"
    const uint32_t version = 1;
    out.write((char*)&magic, sizeof(magic));
    out.write((char*)&version, sizeof(version));

    uint64_t n = mm_chrom_sj.size();
    out.write((char*)&n, sizeof(n));
    for (const auto& [chrom, sjs] : mm_chrom_sj) {
        // write key string
        const uint64_t klen = (uint64_t)chrom.size();
        out.write((char*)&klen, sizeof(klen));
        out.write(chrom.data(), (std::streamsize)klen);

        // write value vector<uint64_t>
        const uint64_t vlen = (uint64_t)sjs.size();
        out.write((char*)&vlen, sizeof(vlen));
        if (vlen) {
            out.write((const char*)sjs.data(), (std::streamsize)(vlen * sizeof(uint64_t)));
        }
    }
}

bool Parser::load_mm_cache_(const fs::path& cache) {
    std::ifstream in(cache, std::ios::binary);
    if (!in) return false;

    uint32_t magic=0, version=0;
    in.read((char*)&magic, sizeof(magic));
    in.read((char*)&version, sizeof(version));
    if (magic != 0x4D4D4348 || version != 1) return false;

    uint64_t n=0;
    in.read((char*)&n, sizeof(n));
    mm_chrom_sj.clear();

    for (uint64_t i=0; i<n; ++i) {
        // read key string
        uint64_t klen=0;
        in.read((char*)&klen, sizeof(klen));
        std::string chrom;
        chrom.resize((size_t)klen);
        if (klen) in.read(chrom.data(), (std::streamsize)klen);

        // read value vector<uint64_t>
        uint64_t vlen=0;
        in.read((char*)&vlen, sizeof(vlen));
        std::vector<uint64_t> sjs;
        sjs.resize((size_t)vlen);
        if (vlen) {
            in.read((char*)sjs.data(), (std::streamsize)(vlen * sizeof(uint64_t)));
        }

        mm_chrom_sj.emplace(std::move(chrom), std::move(sjs));
    }
    return true;
}

//instead of read_mm(filename)
void Parser::read_mm_cached_always(const std::string& filename) {
    fs::path mm_path = filename;
    fs::path cache   = mm_cache_path(mm_path);

    if (fs::exists(cache)) {
        std::cout << "Loading MM cache: " << cache << std::endl;
        if (load_mm_cache_(cache)) return;           // done
        std::cout << "Cache corrupt; rebuilding…" << std::endl;
    }

    // First run (no cache) or cache unreadable: parse once, then save
    std::cout << "Parsing MM and creating cache…" << std::endl;
    read_mm(filename);               // fills mm_chrom_sj
    save_mm_cache_(cache);
    std::cout << "Cached at: " << cache << std::endl;
}

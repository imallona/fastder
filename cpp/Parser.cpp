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
        // calculate total number of reads that map to this bp interval
        // TODO think about int -> unsigned int type safety
        row.length = row.end - row.start;
        // end is not inclusive, since row1.end == row2.start of the next row
        row.total_reads = row.length * row.coverage; //if start = 22, end = 25, coverage = 3 --> (25 - 22) * 3 = 3 * 3 = 9
        library_size += row.total_reads;
        //row.print();
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
            //std::cout << line << std::endl;
            continue;
        }
        SJRow row = SJRow();
        iss >> row;

        //std::cout << row << std::endl;

        rr_all_sj.push_back(row);
    }


}


// pass dictionary all_mm by reference and fill it with keys (sj_id) and values (cumulative count of this sj across samples)
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
            unsigned int sj_id, mm_id, count;

            if (!(iss >> sj_id >> mm_id >> count)){
                //std::cout << "malformed line in MM file: " << line << std::endl;
                continue;
            }

            // OLD: mm_by_samples[sample_id].push_back(std::make_pair(sj_id, count));
            // NEW: cumulative counts of a sj_id across all samples in the input

            // find the rail_id based on the mm_id --> only add counts if the mm_id is part of the samples

            auto it = std::find_if(rail_id_to_mm_id.begin(), rail_id_to_mm_id.end(), [&] (const auto& p)
            {
                return p.second == mm_id;
            });

            // add count if the mm was found
            if (it != rail_id_to_mm_id.end())
            {
                //std::cout << it->first << " : " << it->second << std::endl;
                mm_sj_counts[sj_id] += count; // this creates the binding if it doesn't exist yet, initializes it to 0 and then increases it by count

            }
        }


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
            //convert to integer
            int rail_id = std::stoi(rail_id_str);
            rail_id_to_ext_id.push_back(std::make_pair(rail_id, sample_id));
        }


    }

    // sorting is n log n and finding the position + inserting can be n*n, so better to push_back and then sort


}

// creates a map of rail_id to mm id in rail_id_to_mm_id
// bedgraph_files contains the file name
void Parser::fill_up(std::vector<std::string> bedgraph_files)
{
    //fill up rail_id_to_mm_id
    for (auto& bedgraph_file : bedgraph_files)
    {
        // add the sample and its mm index (= the rank of the rail id across the study) to rail_id_to_mm
        // [&] references all necessary variables i.e. the required context, here it's filename
        auto it = std::find_if(rail_id_to_ext_id.begin(), rail_id_to_ext_id.end(), [&](const auto& sample)
        {
            // the external id is part of the filename for all three sources GTEX, TCGA and SRA
            return bedgraph_file.find(sample.second) != std::string::npos;
        });

        unsigned int mm_id = std::distance(rail_id_to_ext_id.begin(), it) + 1; // std::distance counts the steps between two iterators --> mm_id is 1 too small
        unsigned int rail_id = it->first;

        rail_id_to_mm_id.push_back(std::make_pair(rail_id, mm_id));

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

        if (filename.find("BigWig_list") != std::string::npos && filename.find(".csv") != std::string::npos) //TODO I checked some filenames of the URL csv files manually and they all contain the substring BigWig_list, so I hope that this is a general rule
        {
            std::cout << "BigWig URL list" << std::endl;
            read_url_csv(filename);
            contains_ids = true;
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

    std::cout << "nr of samples in this study:  " << rail_id_to_ext_id.size() << std::endl;

    // fill up rail_id_to_mm_id mapping for all rail_ids in the dataset
    fill_up(bedgraph_files);
    std::cout << "nr of samples provided by user: " << rail_id_to_mm_id.size() << std::endl;

    // now parse all other files
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::string filename = entry.path().string();

        //std::cout << filename << std::endl;
        // read RR file
        if (filename.find("ALL.RR") != std::string::npos) {
            std::cout << "RR file" << std::endl;
            read_rr(filename);

        }
        // don't read in cached MM files as regular MM files!
        else if (entry.path().extension().string() == ".MM" && filename.find("ALL.MM") != std::string::npos && filename.find("mmcache") == std::string::npos ) {
            std::cout << "MM file"<< std::endl;
            // TODO change back read_mm(filename);
            read_mm_cached_always(filename);

        }

        else if (filename.find(".bedGraph") != std::string::npos) {
            std::cout << "Bedgraph file"<< std::endl;

            uint64_t library_size = 0; // ensure that the number is large enough
            std::vector<BedGraphRow> sample_bedgraph = read_bedgraph(filename, library_size);

            //normalize to CPM

            for (BedGraphRow& row : sample_bedgraph)
            {
                row.normalize(library_size);
                //TODO here compute_per_base_coverage(row)
            }

            // add to matrix of all bedgraphs per sample
            all_bedgraphs.push_back(sample_bedgraph);

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

        // create per base coverage across all samples (sample-agnostic)

        //get_per_base_coverages();

    }

    std::cout << "nr of splice junctions across all samples in user input: " << mm_sj_counts.size() << std::endl;
    // for (auto& it : mm_sj_counts)
    // {
    //     std::cout << it.first << " : " << it.second << std::endl;
    // }
    std::cout << "FINISHED PARSING" << std::endl;

}


// CACHE STUFF

namespace fs = std::filesystem;

// Pick where to store the cache (next to the source, simple)
static fs::path mm_cache_path(const fs::path& src) {
    //return mm_path.parent_path() / (mm_path.filename().string() + ".mmcache");
    if (src.extension() == ".mmcache") return src;      // already a cache
    fs::path cache = src;
    cache += ".mmcache";                                // keep original name + add suffix
    // Alternatively: cache.replace_extension(".mmcache"); // if you want to replace .MM
    return cache;
}

// --- serialize your parsed result (mm_sj_counts) ---
void Parser::save_mm_cache_(const fs::path& cache) const {
    std::ofstream out(cache, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open cache for write: " + cache.string());

    // small header: magic + version (helps future-proofing)
    const uint32_t magic = 0x4D4D4348; // "MMCH"
    const uint32_t version = 1;
    out.write((char*)&magic, sizeof(magic));
    out.write((char*)&version, sizeof(version));

    uint64_t n = mm_sj_counts.size();
    out.write((char*)&n, sizeof(n));
    for (const auto& [sj_id, count] : mm_sj_counts) {
        out.write((char*)&sj_id,  sizeof(sj_id));   // unsigned int
        out.write((char*)&count,  sizeof(count));   // unsigned int
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
    mm_sj_counts.clear();
    mm_sj_counts.reserve(n);
    for (uint64_t i=0; i<n; ++i) {
        unsigned int key, val;
        in.read((char*)&key, sizeof(key));
        in.read((char*)&val, sizeof(val));
        mm_sj_counts[key] = val;
    }
    return true;
}

// Call this instead of read_mm(filename)
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
    read_mm(filename);               // fills mm_sj_counts
    save_mm_cache_(cache);
    std::cout << "Cached at: " << cache << std::endl;
}

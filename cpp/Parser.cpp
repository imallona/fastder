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
#include <thread>
#include <future>
#include <cstdint> // for library size which can be too large for unsigned int

#include "Parser.h"

#include <filesystem>

// constructor
Parser::Parser(std::string path_, std::vector<std::string> chromosomes_) {
    path = path_;
    // default: use all chromosomes
    if (chromosomes_.empty())
    {
        std::cout << "[INFO] No chromosomes provided!" << std::endl;
        chromosomes = permitted_chromosomes;
    }

    else
    {
        for (auto chr : chromosomes_)
        {
            // add chr to list of whitelisted chromosomes
            if (std::find(permitted_chromosomes.begin(), permitted_chromosomes.end(), chr) != permitted_chromosomes.end())
            {
                chromosomes.push_back(chr);
            }
        }
    }

}

// function to prevent using splice junctions for chromosomes that weren't parsed
bool Parser::chr_permitted(const std::string& chr) const
{
    for (auto& provided_chr : this->chromosomes)
    {
        if (chr == provided_chr)
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
    //std::cout << "[FILE] " << filename << std::endl;
    //read in file from path
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "[ERROR] could not open .bedgraph file " << filename << std::endl;
    }
    std::string line;
    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);
        //std::cout << line  << std::endl;
        BedGraphRow row = BedGraphRow();
        iss >> row.chrom >> row.start >> row.end >> row.coverage;
        // check if the row is part of the chromosome list passed by the user
        if (std::find(chromosomes.begin(), chromosomes.end(), row.chrom) != chromosomes.end()){
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
    }
    std::cout << "[INFO] Total number of splice junctions: " << rr_all_sj.size() << std::endl;
    //assert(rr_all_sj.size() == 9484210);

}


// create dictionary mm_sj_counts with keys (sj_id) and values (cumulative count of this sj across samples)
// IMPORTANT: the RR file is NOT sorted by chromosomes!
void Parser::read_mm(std::string filename) {
        //read in file from path
        std::ifstream file(filename);
        // max index is 2931 (= nr of samples)
        // min index is 0
        // mm_by_samples.size() = 2931

        if (!file.is_open())
        {
            std::cerr << "[ERROR] could not open file " << filename << std::endl;
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
                //std::cout << nr_of_sj << ", " << rr_all_sj.size() << std::endl;
                if (nr_of_sj != rr_all_sj.size()) {
                    std::cerr << "[ERROR] RR File and number of splice junctions are not equal! Quitting...";
                    return;
                }
                seen_header = true;
                continue;
            }

            // skip invalid lines
            uint64_t sj_id;
            unsigned int mm_id, count;

            if (!(iss >> sj_id >> mm_id >> count)){
                std::cout << "[ERROR] Malformed line in MM file: " << line << std::endl;
                std::cout << line << std::endl;
                continue;
            }

            // OLD: mm_by_samples[sample_id].push_back(std::make_pair(sj_id, count));
            // NEW: cumulative counts of a sj_id across all samples in the input

            // find the rail_id based on the mm_id --> only add sj_id if the mm_id is part of the samples
            auto it = std::find_if(rail_id_to_mm_sample_id.begin(), rail_id_to_mm_sample_id.end(), [&] (const auto& p)
            {
                return p.second == mm_id;
            });

            // add count if the mm was found and if chr is in bedgraph
            //std::cout << rr_all_sj[sj_id].chrom << " for splice junction " << sj_id << std::endl;

            if (it != rail_id_to_mm_sample_id.end() && this->chr_permitted(rr_all_sj[sj_id - 1].chrom)) // rail_id_to_mm_id has <rail_id, mm_id> mapping
            {
                // store vector of sj_ids for each chromosome
                mm_chrom_sj[rr_all_sj[sj_id - 1].chrom].push_back(sj_id); // this creates the binding if it doesn't exist yet, initializes it to 0 and then increases it by count
                //std::cout << mm_chrom_sj[rr_all_sj[sj_id - 1].chrom].back() << std::endl;
            }
            //
        }
        std::cout << "[INFO] MM file contains " << count_lines << " lines"<< std::endl;
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
        std::cerr << "[ERROR] Could not open file " << filename << std::endl;
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
	    //std::cout << bedgraph_file << std::endl;
	    //std::cout << "[" << rail_id_to_ext_id.begin()->second << "]" << std::endl;
        // add the sample and its mm_id (= the rank of the rail id across the study, so all files in total) to rail_id_to_mm
        // [&] references all necessary variables i.e. the required context, here it's filename
        auto it = std::find_if(rail_id_to_ext_id.begin(), rail_id_to_ext_id.end(), [&](auto& sample)
        {
            // search for the external_id in rail_id_to_ext_id and then obtain the rail_id
            // the external id is part of the filename for all three sources GTEX, TCGA and SRA
            //std::cout <<" target " << sample.second << ", size =" << sample.second.size() << std::endl;
            //std::cout <<" bedgraph file " << bedgraph_file << std::endl;
	        sample.second.erase(std::remove(sample.second.begin(), sample.second.end(), '"'),
	            sample.second.end());
            return bedgraph_file.find(sample.second) != std::string::npos;
        });
        if (it != rail_id_to_ext_id.end())
        {
            unsigned int mm_id = std::distance(rail_id_to_ext_id.begin(), it) + 1; // std::distance counts the steps between two iterators --> mm_id is 1 too small, so add 1
            unsigned int rail_id = it->first;

            rail_id_to_mm_sample_id.push_back(std::make_pair(rail_id, mm_id));
        }
        else
        {
            std::cout << "[ERROR] File " << bedgraph_file << " has no rail_id! " << std::endl;
        }
    }
}

void Parser::read_all_bedgraphs(std::vector<std::string> bedgraph_files, unsigned int nof_threads) {
    std::cout << "[INFO] Using " << nof_threads << " threads" << std::endl;
    // reserve space
    all_bedgraphs.resize(bedgraph_files.size());
    all_per_base_coverages.resize(bedgraph_files.size());

    // storage for threads
    std::vector<std::thread> threads;
    threads.reserve(nof_threads);
    // mutex to write to file
    static std::mutex mutex;

    // threads sanity check
    if (nof_threads > bedgraph_files.size()) {
        std::cerr << "[ERROR] Too many threads! Quitting...";
        return;
    }

    // atomic number for index --> never shared, so each index is used exactly once
    // sequence of samples within all_bedgraphs is irrelevant
    std::atomic_int next_index{0};

    for (unsigned int t = 0; t < nof_threads; ++t) {
        threads.emplace_back([this, &bedgraph_files, &next_index]() {
            // infinite loop to ensure that each thread takes the next bedgraph in the queue when it's done
            while (true) {
                unsigned int i = next_index++; //passes index, then does post-increment!
                if (i >= bedgraph_files.size()) break;

                // scope to ensure print statement is not shuffled from concurrency
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    std::cout << "[FILE] Processing BedGraph File " << bedgraph_files.at(i) << std::endl;
                }
                uint64_t library_size = 0; // ensure that the integer type is large enough
                std::vector<BedGraphRow> sample_bedgraph = read_bedgraph(bedgraph_files.at(i), library_size);
                std::unordered_map<std::string, std::vector<double>> per_base_coverage;
                //std::cout << "[INFO] Library size: " << library_size<< std::endl;

                //normalize to CPM and expand rows to per base coverage (also normalized)
                for (BedGraphRow& row : sample_bedgraph)
                {
                    row.normalize(library_size);
                    compute_per_base_coverage(row, per_base_coverage);
                    // row.print();
                }

                // add to matrix of all bedgraphs per sample
                all_bedgraphs[i] = std::move(sample_bedgraph);
                all_per_base_coverages[i] = std::move(per_base_coverage);
            }
        });
    }

    for (auto& thr: threads) {
        thr.join();
    }
}

// attempt to parse all files in path (not recursive!)
void Parser::search_directory() {

    bool contains_ids = false;

    // first check for the external_id to rail_id mapping CSV file
    std::vector<std::string> bedgraph_files;
    std::string mm_file;
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::string filename = entry.path().string();
        // create rail_id_to_ext_id
        if (filename.find("BigWig_list") != std::string::npos && filename.find(".csv") != std::string::npos) //TODO I checked some filenames of the URL csv files manually and they all contain the substring BigWig_list, so I hope that this is a general rule
        {
            std::cout << "[INPUT] BigWig URL list " << filename << std::endl;
            read_url_csv(filename);
            contains_ids = true;
        }
        // read RR file
        else if (filename.find("ALL.RR") != std::string::npos) {
            std::cout << "[INPUT] RR file" << filename << std::endl;
            read_rr(filename);

        }

        // collect all bedgraph files to later fill up rail_id_to_mm_id
        else if (filename.find(".bedGraph") != std::string::npos)
        {
            std::cout << "[INPUT] Bedgraph file "<< filename << std::endl;
            bedgraph_files.push_back(filename);
        }

        else if (entry.path().extension().string() == ".MM" && filename.find("ALL.MM") != std::string::npos && filename.find("mmcache") == std::string::npos) {
            std::cout << "[INPUT] MM file " << filename << std::endl;
            mm_file = filename;
        }
        else {
            std::cout << "[INFO] Unknown file category: " << filename  << std::endl;
        }
    }

    // program cannot run with missing BigWig URL list
    if (!contains_ids || mm_file.empty() || bedgraph_files.empty())
    {
        std::cerr << "[ERROR] Missing input file! Exiting...";
        return;
    }

    // sort rail_id_to_ext_id by rail_id to obtain the index used in the MM file
    std::sort(rail_id_to_ext_id.begin(), rail_id_to_ext_id.end(), [](const auto& a, const auto& b)
    {
        return a.first < b.first;
    });

    std::cout << "[INFO] The study contains " << rail_id_to_ext_id.size() << " samples. " << std::endl;

    // fill up rail_id_to_mm_id mapping for all rail_ids provided by the user
    fill_up(bedgraph_files);
    std::cout << "[INFO] User provided " << rail_id_to_mm_sample_id.size() << " bedgraph files." << std::endl;

    unsigned int max_threads = std::max(int(std::thread::hardware_concurrency()), 4); //require at least 4 cores
    unsigned int nof_samples =  rail_id_to_mm_sample_id.size();
    unsigned int nof_threads = std::min(max_threads, nof_samples);

    // launch separate thread to parse MM file
    std::thread mm_thread(&Parser::read_mm, this, mm_file);

    // parse all bedgraph files concurrently
    read_all_bedgraphs(bedgraph_files, nof_threads);

    // stop MM thread
    mm_thread.join();
    std::cout << "[INFO] Finished parsing all files." << std::endl;

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

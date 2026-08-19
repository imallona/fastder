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
#include <atomic>
#include <mutex>
#include <charconv>
#include <cstdint> // for library size which can be too large for unsigned int
#include <cstdlib>

#include "Parser.h"

// constructor
Parser::Parser(std::string path_, std::vector<std::string> chromosomes_, int cores_) {
    path = path_;
    user_cores = cores_;
    std::cout << "[INFO] fastder will use up to " << cores_ << " cores. To change the maximum number of cores, provide a different value with the --cores flag." << std::endl;

    if (chromosomes_.empty())
    {
        std::cout << "[INFO] User specified no chromosomes. fastder uses all chromosomes by default.\n" << std::endl;
        chromosomes_vec.assign(permitted_chromosomes.begin(), permitted_chromosomes.end());
        chromosomes_set = permitted_chromosomes;

    }

    else
    {
        for (auto chr : chromosomes_)
        {
            // add chr to list of whitelisted chromosomes
            if (permitted_chromosomes.contains(chr))
            {
                chromosomes_vec.emplace_back(chr);
                chromosomes_set.insert(chr);
            }
        }
    }
    // if only invalid chromosomes were provided
    if (chromosomes_vec.empty() || chromosomes_set.empty())
    {
        std::cout << "[INFO] User specified no valid chromosomes! fastder will use all chromosomes." << std::endl;
        chromosomes_vec.assign(permitted_chromosomes.begin(), permitted_chromosomes.end());
        chromosomes_set = permitted_chromosomes;
    }
}

// compute_per_base_coverage was deleted. The dense per-base double vector it
// produced was the largest single resident structure in fastder (47 Mb x 8 B
// per chr21 sample, 24 GB per full-hg38 sample). Averager now consumes the
// sparse interval form in all_bedgraphs directly.

// parse relevant chromosomes of a bedgraph file
std::vector<BedGraphRow> Parser::read_bedgraph(const std::string& filename, uint64_t& library_size) const
{
    std::vector<BedGraphRow> bedgraph; // stores the full bedgraph of one sample, organized by rows (bins) with the same coverage

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
        if (line.empty()) return bedgraph;

        BedGraphRow row = BedGraphRow();
        // read in line by line with std::from_chars()
        const char* p = line.data();
        const char* end = p + line.size(); // end of the line

        while (p < end &&  (*p == ' ' || *p == '\t')) p++; // skip ws or tab at the beginning
        const char* start = p;
        while (p < end && (*p != ' ' && *p != '\t')) ++p;
        std::string_view chrom(start, p - start);

        // skip whitespace and tab
        while (p < end &&  (*p == ' ' || *p == '\t')) p++;

        std::from_chars_result res1 = std::from_chars(p, end, row.start);
        if (res1.ec != std::errc{})
        {
            std::cerr << "[ERROR] Malformed start position in Bedgraph file: " << line << std::endl;
            continue;
        }
        p = res1.ptr;

        // skip whitespace and tab
        while (p < end &&  (*p == ' ' || *p == '\t')) p++;

        std::from_chars_result res2 = std::from_chars(p, end, row.end);
        if (res2.ec != std::errc{})
        {
            std::cerr << "[ERROR] Malformed end position in Bedgraph file: " << line << std::endl;
            continue;
        }
        p = res2.ptr;

        // skip whitespace and tab
        while (p < end &&  (*p == ' ' || *p == '\t')) p++;

        int coverage_as_int; //std::from_chars only exists for integers and initially, column 4 is of type integer in Bedgraphs
        std::from_chars_result res3 = std::from_chars(p, end, coverage_as_int);
        if (res3.ec != std::errc{})
        {
            std::cerr << "[ERROR] Malformed coverage entry in Bedgraph file: " << line << std::endl;
            continue;
        }
        row.coverage = static_cast<float>(coverage_as_int);

        // calculate total number of reads that map to this bp interval
        row.length = row.end - row.start;
        // end is not inclusive, since row1.end == row2.start of the next row
        row.total_reads = static_cast<uint64_t>(row.length * row.coverage); //if start = 22, end = 25, coverage = 3 --> (25 - 22) * 3 = 3 * 3 = 9

        // Library size covers the whole file, not only the chromosomes in
        // --chr. Summing over the selection instead would make the same
        // --min-coverage a different absolute cutoff in every run. Artificial
        // chromosomes (chrUn_GL000218v1 and the like) are part of the library
        // too, so they are counted here and dropped below, where canonical
        // names are the ones 4 or 5 characters long.
        library_size += row.total_reads;

        if (chrom.size() != 4 && chrom.size() != 5) {
            continue;
        }
        row.chrom = std::string(chrom);

        // check if the row is part of the chromosome list passed by the user
        if (chromosomes_set.contains(row.chrom)){
            bedgraph.emplace_back(row);
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
        std::cerr << "[ERROR] Could not open .rr file " << filename << std::endl;
    }
    std::string line;

    // Two-pass parse: keep only RR rows whose chromosome is in
    // chromosomes_set. The on-disk MM references rows by their 1-based RR
    // row number, so we record a remap from old sj_id to the new compact
    // index in rr_all_sj. Dropped rows get 0 in the remap.
    // For full hg38 with --chr chr21, this drops rr_all_sj from about 9.5M
    // rows to about 3.5k rows. The remap itself costs 4 bytes per RR row.
    rr_all_sj.clear();
    sj_id_remap.clear();
    rr_total_rows = 0;

    while (std::getline(file, line))
    {
        // skip invalid lines and headers (which contain the string "chromosome")
        if (line.empty() || line.find("chromosome") != std::string::npos) {
            continue;
        }

        SJRow row;
        std::istringstream iss(line);
        iss >> row;

        ++rr_total_rows;

        if (chromosomes_set.contains(row.chrom)) {
            rr_all_sj.emplace_back(std::move(row));
            // 1-based index into rr_all_sj
            sj_id_remap.push_back(static_cast<uint32_t>(rr_all_sj.size()));
        } else {
            sj_id_remap.push_back(0);  // sentinel: row not retained
        }
    }
    std::cout << "[INFO] RR rows total: " << rr_total_rows
              << " | retained for analysed chromosomes: " << rr_all_sj.size()
              << std::endl;
}



// read MM file
// IMPORTANT: the RR file is not sorted by chromosomes!
void Parser::read_mm(std::string filename) {
        // Reset MM-derived state in case the same Parser instance reads more
        // than one MM file. Junction ids would otherwise accumulate across
        // calls and downstream stitching would see duplicates.
        mm_chrom_sj.clear();

        //read in file from path
        std::ifstream file(filename);
        //
        // static thread_local std::vector<char> buf(1 << 20); // set buffer to 1mb;
        // file.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
        // max index is 2931 (= nr of samples)
        // min index is 0
        // mm_by_samples.size() = 2931

        if (!file.is_open())
        {
            std::cerr << "[ERROR] could not open MM file " << filename << std::endl;
        }
        // Read support per junction, summed over the loaded samples. Only
        // filled when a threshold is set: at 8 bytes per junction this costs
        // about 76 MB on the 9.5M junctions of a full hg38 RR file.
        std::vector<uint64_t> junction_reads;
        if (min_junction_reads > 0) junction_reads.assign(rr_all_sj.size(), 0);

        std::string line;
        bool seen_header = false;
        uint64_t nr_of_sj = 0;
        uint64_t sj_occ_in_samples = 0;
        unsigned int nr_of_samples = 1;
        //auto sj_id_prev = 0;
        uint64_t count_lines = 0;
        while (std::getline(file, line))
        {
            ++count_lines;
            // read in line by line
            if (line.empty()) return;
            if (line[0] == '%') continue;
            if (!seen_header) {
                std::istringstream iss(line);
                // header: 9484210	2931	699368828, actual #lines = 699368831
                iss >> nr_of_sj >> nr_of_samples >> sj_occ_in_samples;
                // Compare against the *total* number of rows seen in RR, not
                // rr_all_sj.size(): rr_all_sj is now a chr-filtered subset.
                if (nr_of_sj != rr_total_rows) {
                    std::cerr << "[ERROR] MM header sj count (" << nr_of_sj
                              << ") does not match RR row count (" << rr_total_rows
                              << "). Quitting..." << std::endl;
                    return;
                }
                seen_header = true;
                continue;
            }

            uint64_t sj_id = 0;
            unsigned int mm_id = 0;
            unsigned int count = 0;
            // use std::from_chars for faster, manual parsing
            const char* p = line.data();
            const char* end = p + line.size(); // end of the line


            std::from_chars_result res1 = std::from_chars(p, end, sj_id);
            if (res1.ec != std::errc{})
            {
                std::cerr << "[ERROR] Malformed line in MM file: " << line << std::endl;
                continue;
            }
            p = res1.ptr;

            // skip whitespace and tab
            while (p < end &&  (*p == ' ' || *p == '\t')) p++;

            std::from_chars_result res2 = std::from_chars(p, end, mm_id);
            if (res2.ec != std::errc{})
            {
                std::cerr << "[ERROR] Malformed line in MM file: " << line << std::endl;
                continue;
            }

            // Look up the new (post-filter) sj_id via the remap built in
            // read_rr. A 0 means this junction was on a chromosome we don't
            // analyse and was dropped, so skip without ever touching rr_all_sj.
            if (sj_id == 0 || sj_id - 1 >= sj_id_remap.size()) continue;
            const uint32_t new_sj_id = sj_id_remap[sj_id - 1];
            if (new_sj_id == 0) continue;

            if (mm_ids.contains(mm_id))
            {
                if (min_junction_reads > 0)
                {
                    // Third column of the MatrixMarket triplet: reads
                    // supporting this junction in this sample. res2 left p
                    // before mm_id, so advance past it first.
                    p = res2.ptr;
                    while (p < end && (*p == ' ' || *p == '\t')) p++;
                    unsigned int count = 0;
                    if (std::from_chars(p, end, count).ec == std::errc{})
                    {
                        junction_reads[new_sj_id - 1] += count;
                    }
                }
                // mm_chrom_sj stores *new* sj_ids (indexes into rr_all_sj).
                mm_chrom_sj[rr_all_sj[new_sj_id - 1].chrom].emplace_back(new_sj_id);
            }
        }

        if (min_junction_reads > 0)
        {
            uint64_t dropped = 0;
            for (auto& [chrom, ids] : mm_chrom_sj)
            {
                auto kept = std::remove_if(ids.begin(), ids.end(),
                    [&](uint32_t id) { return junction_reads[id - 1] < min_junction_reads; });
                dropped += static_cast<uint64_t>(std::distance(kept, ids.end()));
                ids.erase(kept, ids.end());
            }
            std::cout << "[INFO] Dropped " << dropped
                      << " junction entries below --min-junction-reads "
                      << min_junction_reads << std::endl;
        }
        //std::cout << "[INFO] MM file contains " << count_lines << " lines"<< std::endl;
        //assert(sj_occ_in_samples <= count_lines);
    }

// parse bigwig URL list csv file
void Parser::read_url_csv(std::string filename)
{
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

        // invalid or header line --> the first line is always skipped
        if (line.empty() || line == "rail_id,external_id,study,BigWigURL" || first_line)
        {
            first_line = false;
            continue;
        }
        std::string rail_id_str, sample_id;
        // only read in the first two tab-separated entries, ignore the rest!
        if (std::getline(iss, rail_id_str, ',') && std::getline(iss, sample_id, ','))
        {
            //convert to integer
            int rail_id = std::stoi(rail_id_str);
            rail_id_to_ext_id.emplace_back(std::make_pair(rail_id, sample_id));
        }
    }
    // sorting is n log n and finding the position + inserting can be n*n, so better to emplace_back and then sort
}

// creates a map of rail_id to mm_id in rail_id_to_mm_id
// bedgraph_files contains the file names of all samples
void Parser::fill_up(std::vector<std::string> bedgraph_files)
{
    //fill up rail_id_to_mm_id
    for (auto& bedgraph_file : bedgraph_files)
    {
        // add the sample and its mm_id (= the rank of the rail id across the study, so all files in total) to rail_id_to_mm
        // [&] references all necessary variables i.e. the required context (here: filename)
        auto it = std::find_if(rail_id_to_ext_id.begin(), rail_id_to_ext_id.end(), [&](auto& sample)
        {
            // search for the external_id in rail_id_to_ext_id and then obtain the rail_id
            // the external id is part of the filename for all three sources GTEX, TCGA and SRA
	        sample.second.erase(std::remove(sample.second.begin(), sample.second.end(), '"'),
	            sample.second.end());
            return bedgraph_file.find(sample.second) != std::string::npos;
        });
        if (it != rail_id_to_ext_id.end())
        {
            unsigned int mm_id = std::distance(rail_id_to_ext_id.begin(), it) + 1; // std::distance counts the steps between two iterators --> mm_id is 1 too small, so add 1
            mm_ids.insert(mm_id);
        }
        else
        {
            std::cerr << "[ERROR] File " << bedgraph_file << " has no rail_id! Check that the external_id is contained in the sample file name. " << std::endl;
        }
    }
}

void Parser::read_all_bedgraphs(std::vector<std::string> bedgraph_files, unsigned int nof_threads) {
    std::cout << "[INFO] fastder is using " << nof_threads + 1 << " threads for parsing." << std::endl;
    // reserve space
    all_bedgraphs.resize(bedgraph_files.size());

    // storage for threads
    std::vector<std::thread> threads;
    threads.reserve(nof_threads);
    // mutex to write to file
    static std::mutex mutex;

    // atomic number for index --> never shared, so each index is used exactly once
    // sequence of samples within all_bedgraphs is irrelevant
    std::atomic_int next_index{0};

    for (unsigned int t = 0; t < nof_threads; ++t) {
        threads.emplace_back([this, &bedgraph_files, &next_index]() {
            // infinite loop to ensure that each thread takes the next bedgraph in the queue when it's done
            while (true) {
                unsigned int i = next_index++; //passes index, then does post-increment!
                if (i >= bedgraph_files.size()) break;

                const std::string& filename = bedgraph_files.at(i);
                // mutex to ensure print statement is not shuffled from concurrency
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    std::cout << "[FILE] Processing coverage file " << filename << std::endl;
                }
                uint64_t library_size = 0; // ensure that the integer type is large enough
                std::vector<BedGraphRow> sample_bedgraph;

                // pick the right reader by file extension. BigWig parsing is
                // gated on libBigWig at compile time; without it a .bw input
                // is a hard error so we don't silently dilute the mean
                // coverage by counting an empty sample toward total_samples.
                if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".bw")
                {
#ifdef FASTDER_USE_LIBBIGWIG
                    sample_bedgraph = read_bigwig(filename, library_size);
#else
                    std::cerr << "[ERROR] BigWig input " << filename
                              << " requires fastder built with -DFASTDER_USE_LIBBIGWIG=ON. "
                              << "Reconfigure or convert the file to BedGraph." << std::endl;
                    std::exit(1);
#endif
                }
                else
                {
                    sample_bedgraph = read_bedgraph(filename, library_size);
                }

                // normalize each interval to CPM. Per-base expansion is no
                // longer performed; Averager consumes the sparse intervals.
                for (BedGraphRow& row : sample_bedgraph)
                {
                    row.normalize(library_size);
                }

                // add the sample's intervals to the per-sample matrix
                {
                    std::lock_guard lock(mutex);
                    all_bedgraphs[i] = std::move(sample_bedgraph);
                }

            }
        });
    }

    for (auto& thr: threads) {
        thr.join();
    }
}


// libBigWig integration. The body is gated on FASTDER_USE_LIBBIGWIG so the
// default build is hermetic. When the option is on, libBigWig is FetchContent'd
// from upstream and its headers are on the include path.
#ifdef FASTDER_USE_LIBBIGWIG
extern "C" {
#include <bigWig.h>
}
#endif

std::vector<BedGraphRow> Parser::read_bigwig(const std::string& filename, uint64_t& library_size,
                                             char strand) const
{
    std::vector<BedGraphRow> intervals;
#ifdef FASTDER_USE_LIBBIGWIG
    // bwInit allocates a process-wide read buffer. Call it once across all
    // threads. bwCleanup is left to process exit; libBigWig's cleanup function
    // is not safe to call while other readers may still be active.
    static std::once_flag bw_init_flag;
    std::call_once(bw_init_flag, []() { bwInit(1 << 17); });

    bigWigFile_t* fp = bwOpen(const_cast<char*>(filename.c_str()), nullptr, "r");
    if (!fp)
    {
        std::cerr << "[ERROR] Could not open BigWig " << filename << std::endl;
        return intervals;
    }

    // Whole-file library size, as in read_bedgraph. The total summary holds
    // this sum already, so no extra intervals are read. A file written without
    // a summary block needs the full walk instead.
    if (fp->hdr->summaryOffset != 0)
    {
        library_size += static_cast<uint64_t>(fp->hdr->sumData);
    }
    else
    {
        for (int64_t k = 0; k < fp->cl->nKeys; ++k)
        {
            const std::string chrom = fp->cl->chrom[k];
            bwOverlappingIntervals_t* o = bwGetOverlappingIntervals(
                fp, const_cast<char*>(chrom.c_str()), 0, fp->cl->len[k]);
            if (!o) continue;
            for (uint32_t i = 0; i < o->l; ++i)
            {
                library_size += static_cast<uint64_t>(
                    static_cast<double>(o->end[i] - o->start[i]) * o->value[i]);
            }
            bwDestroyOverlappingIntervals(o);
        }
    }

    // Iterate the BigWig's chromosomes. Skip any that the user did not
    // request via --chr (chromosomes_set), matching read_bedgraph's filter.
    for (int64_t k = 0; k < fp->cl->nKeys; ++k)
    {
        const std::string chrom = fp->cl->chrom[k];
        const uint32_t chrom_len = fp->cl->len[k];
        if (!chromosomes_set.contains(chrom)) continue;

        // bwGetOverlappingIntervals returns the BigWig's intrinsic intervals
        // that overlap [start, end). Calling it for the entire chromosome
        // yields the file's stored intervals on that chromosome.
        bwOverlappingIntervals_t* o = bwGetOverlappingIntervals(
            fp, const_cast<char*>(chrom.c_str()), 0, chrom_len);
        if (!o) continue;

        intervals.reserve(intervals.size() + o->l);
        for (uint32_t i = 0; i < o->l; ++i)
        {
            BedGraphRow row(chrom,
                            static_cast<uint64_t>(o->start[i]),
                            static_cast<uint64_t>(o->end[i]),
                            static_cast<double>(o->value[i]),
                            strand);
            row.length = static_cast<unsigned int>(row.end - row.start);
            row.total_reads = static_cast<uint64_t>(row.length * row.coverage);
            intervals.emplace_back(std::move(row));
        }
        bwDestroyOverlappingIntervals(o);
    }

    bwClose(fp);
#else
    (void)filename; (void)library_size; (void)strand;
    std::cerr << "[ERROR] read_bigwig was called but fastder was built without "
                 "libBigWig support. Reconfigure with -DFASTDER_USE_LIBBIGWIG=ON "
                 "or feed BedGraph (.bedGraph) input." << std::endl;
#endif
    return intervals;
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
            std::cout << "[FILE] Processing Metadata CSV File " << filename << std::endl;
            read_url_csv(filename);
            contains_ids = true;
        }
        // read RR file
        else if (filename.find("ALL.RR") != std::string::npos) {
            std::cout << "[FILE] Processing RR File " << filename << std::endl;
            read_rr(filename);

        }

        // collect all coverage files (BedGraph or BigWig) to later fill up
        // rail_id_to_mm_id. The variable name keeps "bedgraph_files" for
        // continuity but the list is just paths; read_all_bedgraphs picks the
        // right reader by extension. .bw files require a libBigWig-enabled build.
        else if (filename.find(".bedGraph") != std::string::npos ||
                 (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".bw"))
        {
            bedgraph_files.emplace_back(filename);
        }

        else if (entry.path().extension().string() == ".MM" && filename.find("ALL.MM") != std::string::npos && filename.find("mmcache") == std::string::npos) {
            mm_file = filename;
        }
        // else {
        //     std::cout << "[INFO] Unknown file category: " << filename  << std::endl;
        // }
    }

    // program cannot run with missing BigWig URL list
    if (!contains_ids || mm_file.empty() || bedgraph_files.empty())
    {
        std::cerr << "[ERROR] Missing input file! Exiting..." << std::endl;
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
    std::cout << "[INFO] User provided " << mm_ids.size() << " samples." << std::endl;

    unsigned int nof_samples =  mm_ids.size();
    unsigned int nof_threads = std::min(user_cores, nof_samples);

    // launch separate thread to parse MM file
    std::cout << "[FILE] Processing MM File " << mm_file << std::endl;
    std::thread mm_thread(&Parser::read_mm, this, mm_file);

    // parse all bedgraph files concurrently
    read_all_bedgraphs(bedgraph_files, nof_threads);

    // stop MM thread
    mm_thread.join();

}

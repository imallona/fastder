//
// Created by martinalavanya on 24.09.25.
//

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include "BedGraphRow.h"
#include "SJRow.h"
#include "Averager.h"


Averager::Averager(int threads_)
{
    nof_threads = threads_;
}

namespace {

// Bucket each sample's intervals by chromosome so the per-chromosome sweep
// only sees its own data. Within each chromosome the intervals are sorted by
// start position because read_bedgraph / read_bigwig already emit in file
// order, but we re-sort defensively in case a future input source mixes order.
std::unordered_map<std::string, std::vector<std::vector<BedGraphRow>>>
bucket_by_chrom(const std::vector<std::vector<BedGraphRow>>& all_bedgraphs)
{
    std::unordered_map<std::string, std::vector<std::vector<BedGraphRow>>> chrom_samples;
    for (const auto& sample : all_bedgraphs)
    {
        // collect this sample's rows per chromosome, then push into the global map
        std::unordered_map<std::string, std::vector<BedGraphRow>> per_chrom;
        for (const auto& row : sample)
        {
            per_chrom[row.chrom].push_back(row);
        }
        for (auto& [chrom, rows] : per_chrom)
        {
            std::sort(rows.begin(), rows.end(),
                      [](const BedGraphRow& a, const BedGraphRow& b) { return a.start < b.start; });
            chrom_samples[chrom].push_back(std::move(rows));
        }
    }
    return chrom_samples;
}

// Compute the sparse mean for one chromosome. samples[s] is sample s's sorted,
// non-overlapping interval list on this chromosome. total_samples is the
// total sample count across the experiment, including samples that have no
// coverage on this chromosome (they contribute zero to the mean).
//
// Algorithm: collect every distinct start/end position from every sample as a
// breakpoint. Between two consecutive breakpoints no sample's interval starts
// or ends, so each sample's contribution is constant over the segment. We
// advance a per-sample iterator alongside the sweep so the per-segment
// lookup is amortised O(1). Zero-mean segments are dropped from the output.
std::vector<BedGraphRow> mean_for_chrom(const std::string& chrom,
                                        const std::vector<std::vector<BedGraphRow>>& samples,
                                        size_t total_samples)
{
    std::vector<BedGraphRow> result;
    if (samples.empty() || total_samples == 0) return result;

    // gather and dedupe breakpoints
    std::vector<uint64_t> breakpoints;
    for (const auto& sample : samples)
    {
        breakpoints.reserve(breakpoints.size() + sample.size() * 2);
        for (const auto& row : sample)
        {
            breakpoints.push_back(row.start);
            breakpoints.push_back(row.end);
        }
    }
    std::sort(breakpoints.begin(), breakpoints.end());
    breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()), breakpoints.end());
    if (breakpoints.size() < 2) return result;

    // per-sample index into samples[s]: points to the first interval whose
    // end > current_position, or samples[s].size() once exhausted.
    std::vector<size_t> idx(samples.size(), 0);

    for (size_t b = 0; b + 1 < breakpoints.size(); ++b)
    {
        const uint64_t seg_start = breakpoints[b];
        const uint64_t seg_end   = breakpoints[b + 1];

        double sum = 0.0;
        for (size_t s = 0; s < samples.size(); ++s)
        {
            // skip past intervals that lie entirely before this segment
            while (idx[s] < samples[s].size() && samples[s][idx[s]].end <= seg_start)
            {
                ++idx[s];
            }
            // include this sample's coverage only if its current interval covers seg_start
            if (idx[s] < samples[s].size() &&
                samples[s][idx[s]].start <= seg_start &&
                seg_start < samples[s][idx[s]].end)
            {
                sum += samples[s][idx[s]].coverage;
            }
        }
        const double mean = sum / static_cast<double>(total_samples);
        if (mean > 0.0)
        {
            // adjacent identical-mean segments are coalesced below
            result.emplace_back(BedGraphRow(chrom, seg_start, seg_end, mean));
        }
    }

    // coalesce contiguous intervals that share the same mean coverage. This
    // keeps the interval count tight and makes find_ERs's contiguity check cheap.
    if (result.size() < 2) return result;
    std::vector<BedGraphRow> merged;
    merged.reserve(result.size());
    merged.emplace_back(result[0]);
    for (size_t i = 1; i < result.size(); ++i)
    {
        BedGraphRow& last = merged.back();
        if (last.end == result[i].start && last.coverage == result[i].coverage)
        {
            last.end = result[i].end;
            last.length = static_cast<unsigned int>(last.end - last.start);
        }
        else
        {
            merged.emplace_back(result[i]);
        }
    }
    return merged;
}

} // anonymous namespace


// compute mean coverage as sparse intervals across samples
void Averager::compute_mean_coverage(std::vector<std::vector<BedGraphRow>>& all_bedgraphs)
{
    if (all_bedgraphs.empty())
    {
        std::cerr << "[ERROR] No samples were provided to compute_mean_coverage." << std::endl;
        return;
    }

    auto chrom_samples = bucket_by_chrom(all_bedgraphs);
    const size_t total_samples = all_bedgraphs.size();

    // populate the chroms list once, used by find_ERs to parallelise
    chroms.clear();
    chroms.reserve(chrom_samples.size());
    for (auto& [chrom, _] : chrom_samples) chroms.push_back(chrom);

    std::vector<std::thread> threads;
    threads.reserve(nof_threads);
    std::atomic_int next_index{0};

    std::cout << "[INFO] fastder is using " << nof_threads << " threads for computing mean coverage." << std::endl;

    for (unsigned int t = 0; t < static_cast<unsigned int>(nof_threads); ++t)
    {
        threads.emplace_back([&]()
        {
            while (true)
            {
                unsigned int i = next_index++;
                if (i >= chroms.size()) break;
                const std::string& chrom = chroms[i];
                std::vector<BedGraphRow> intervals = mean_for_chrom(chrom, chrom_samples[chrom], total_samples);
                {
                    std::lock_guard<std::mutex> lock(map_mutex);
                    mean_intervals[chrom] = std::move(intervals);
                }
            }
        });
    }
    for (auto& thr : threads) thr.join();
}


// find ERs as maximal contiguous runs of mean intervals above threshold
void Averager::find_ERs(double threshold, int min_length)
{
    if (mean_intervals.empty())
    {
        std::cerr << "[ERROR] mean_intervals is empty." << std::endl;
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(nof_threads);
    std::atomic_int next_index{0};

    for (unsigned int t = 0; t < static_cast<unsigned int>(nof_threads); ++t)
    {
        workers.emplace_back([&, t] {
            while (true)
            {
                unsigned int i = next_index++;
                if (i >= chroms.size()) break;
                const std::string& chrom = chroms[i];
                // Read-only lookup: operator[] would silently insert a default
                // entry under concurrent access from worker threads, which is
                // undefined behavior on a shared map.
                const auto intervals_it = mean_intervals.find(chrom);
                if (intervals_it == mean_intervals.end())
                {
                    std::cerr << "[ERROR] Missing mean_intervals entry for chromosome: "
                              << chrom << std::endl;
                    continue;
                }
                const std::vector<BedGraphRow>& intervals = intervals_it->second;

                std::vector<BedGraphRow> chrom_ers;

                // running ER state
                bool in_er = false;
                uint64_t er_start = 0;
                uint64_t er_end   = 0;
                double weighted_sum = 0.0;

                auto finalize = [&]()
                {
                    if (in_er && (er_end - er_start) > static_cast<uint64_t>(min_length))
                    {
                        const double avg = weighted_sum / static_cast<double>(er_end - er_start);
                        chrom_ers.emplace_back(BedGraphRow(chrom, er_start, er_end, avg));
                    }
                    in_er = false;
                    weighted_sum = 0.0;
                };

                for (const auto& row : intervals)
                {
                    if (row.coverage > threshold)
                    {
                        // a gap or non-contiguous step also breaks the ER, even if
                        // the new interval is above threshold
                        if (in_er && row.start == er_end)
                        {
                            er_end = row.end;
                            weighted_sum += row.coverage * static_cast<double>(row.end - row.start);
                        }
                        else
                        {
                            finalize();
                            er_start = row.start;
                            er_end   = row.end;
                            weighted_sum = row.coverage * static_cast<double>(row.end - row.start);
                            in_er = true;
                        }
                    }
                    else
                    {
                        finalize();
                    }
                }
                finalize();

                {
                    std::lock_guard<std::mutex> lock(map_mutex);
                    expressed_regions[chrom] = std::move(chrom_ers);
                }
            }
        });
    }
    for (auto& w : workers) w.join();
}

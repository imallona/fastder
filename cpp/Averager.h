//
// Created by marti on 08/10/2025.
//
#pragma once

#include "BedGraphRow.h"
#include <vector>
#include <iostream>
#include <unordered_map>
#include <mutex>
#ifndef FASTDER_AVERAGE_H
#define FASTDER_AVERAGE_H

#endif //FASTDER_AVERAGE_H

class Averager {

    public:
        Averager(int threads_);
        // Compute the per-chromosome mean coverage as a sparse interval list
        // by sweep-merging the per-sample interval lists. The previous
        // implementation built a dense per-base vector<double> per chromosome
        // per sample; on full hg38 that grew to ~24 GB per sample. The sweep
        // produces only as many output intervals as the union of breakpoints
        // across samples requires, with zero coverage segments suppressed.
        void compute_mean_coverage(std::vector<std::vector<BedGraphRow>>& all_bedgraphs);
        // Walk the sparse mean_intervals and emit one BedGraphRow per
        // expressed region: a maximal contiguous run of intervals all above
        // threshold. min_length is in nucleotides on the reference.
        void find_ERs(double threshold, int min_length);

        int nof_threads;
        std::vector<std::string> chroms;
        std::mutex map_mutex;
        // Mean coverage per chromosome as a sorted, non-overlapping list of
        // BedGraphRows. Adjacent rows can be contiguous (row[i].end == row[i+1].start);
        // any uncovered region between rows is implicit zero coverage.
        std::unordered_map<std::string, std::vector<BedGraphRow>> mean_intervals;
        // Expressed regions per chromosome, one BedGraphRow per ER. Strand on
        // each ER is inherited from mean_intervals (currently '.' since the
        // BigWig input is unstranded; the SJ-strand stitching that the
        // Integrator does later is what tags the stitched ERs).
        std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;



};
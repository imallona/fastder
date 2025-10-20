//
// Created by marti on 08/10/2025.
//
#pragma once

#include "BedGraphRow.h"
#include <vector>
#include <iostream>
#include <unordered_map>
#ifndef FASTDER_AVERAGE_H
#define FASTDER_AVERAGE_H

#endif //FASTDER_AVERAGE_H

class Averager {

    public:
        void compute_mean_coverage();
        //bool in_interval(double current_avg, double bp_coverage);
        void find_ERs(double threshold, int min_length);
        void get_all_per_base_coverage(const std::vector<std::vector<BedGraphRow>>& all_bedgraphs);
        static void compute_per_base_coverage(const BedGraphRow& row, std::unordered_map<std::string, std::vector<double>>& per_base_coverage);
        void stitch_up();

        // MEMBER VARIABLES

        //matrix consisting of multiple vectors, where each row in a vector is of type BedGraphRow,
        //a BedGraphRow is a bin of nucleotides with the same read count
        // a vector of BedGraphRows corresponds to one sample
        //matrix where each vector contains the normalized count per base from ONE sample

        std::unordered_map<std::string, std::vector<double>> mean_coverage; //key = chromosome, value = BedGraphRow
        std::vector<BedGraphRow> expressed_regions;
        std::vector<std::unordered_map<std::string, std::vector<double>>> all_per_base_coverages;
        // one unordered map per sample with keys = chromosome nr, values = vector of per-base coverage for that chromosome
        // store all the individual sample maps in a vector (since sample identity doesn't matter anymore later on)



};
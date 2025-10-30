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
        void compute_mean_coverage(std::vector<std::unordered_map<std::string, std::vector<double>>>& all_per_base_coverages, const std::vector<std::string>& chromosome_sequence);

        void find_ERs(double threshold, int min_length, std::vector<std::string>& chromosome_sequence);

        // MEMBER VARIABLES

        //matrix consisting of multiple vectors, where each row in a vector is of type BedGraphRow,
        //a BedGraphRow is a bin of nucleotides with the same read count
        // a vector of BedGraphRows corresponds to one sample
        //matrix where each vector contains the normalized count per base from ONE sample

        std::unordered_map<std::string, std::vector<double>> mean_coverage; //key = chromosome, value = vector of per-base mean coverage of this chromosome
        std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions; //key = chromosome, value = vector of BedGraphRows with coverage > threshold
        //std::vector<std::unordered_map<std::string, std::vector<double>>> all_per_base_coverages;
        // one unordered map per sample with keys = chromosome nr, values = vector of per-base coverage for that chromosome
        // store all the individual sample maps in a vector (since sample identity doesn't matter anymore later on)



};
//
// Created by marti on 08/10/2025.
//
#pragma once

#include "BedGraphRow.h"
#include <vector>
#include <iostream>
#ifndef FASTDER_AVERAGE_H
#define FASTDER_AVERAGE_H

#endif //FASTDER_AVERAGE_H

class Average {

    public:
        void compute_per_base_coverage(const BedGraphRow& row, std::vector<double>& per_base_coverage);
        void normalize(std::vector<double>& per_base_coverage, const int library_size);
        std::vector<double> compute_avg_coverage();
        bool in_interval(double current_avg, double bp_coverage);
        std::vector<BedGraphRow> find_ERs(const std::vector<double>& avg_coverage);



        // MEMBER VARIABLES

        //matrix consisting of multiple vectors, where each row in a vector is of type BedGraphRow,
        //a BedGraphRow is a bin of nucleotides with the same read count
        // a vector of BedGraphRows corresponds to one sample
        std::vector<std::vector<BedGraphRow>> all_bedgraphrows;
        //matrix where each vector contains the normalized count per base from ONE sample
        std::vector<std::vector<double>> all_per_base_coverages;
        std::vector<BedGraphRow> results;

};
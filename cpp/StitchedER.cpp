//
// Created by martinalavanya on 20.10.25.
//

#include "StitchedER.h"

// later change to add weight / memory
double StitchedER::get_avg_coverage(){
    double sum = 0;
    unsigned int total_length = 0;
    for (auto er : this->all_coverages){
        sum += er.first * er.second; // weighted sum of avg coverages across ERs
        total_length += er.second;

    }
    return sum / total_length;
}
// function that calculates relative match with a tolerance of +/- 5%
bool StitchedER::is_similar(double val1, double val2){
    double tolerance_bottom = val1 * 0.95;
    double tolerance_top = val1 * 1.05;
    return val2 >= tolerance_bottom && val2 <= tolerance_top;
}


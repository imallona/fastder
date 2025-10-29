//
// Created by martinalavanya on 24.09.25.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <cassert>
#include "BedGraphRow.h"
#include "SJRow.h"
#include "Averager.h"

//compute overall average coverage
void Averager::compute_mean_coverage(std::vector<std::unordered_map<std::string, std::vector<double>>>& all_per_base_coverages)
{

    //iterate
    std::cout << "#samples = " << all_per_base_coverages.size() << ", #chromosomes = "<< all_per_base_coverages[0].size() << std::endl;

    // iterate over chromosomes
    for (auto& pair : all_per_base_coverages[0])
    {
        std::string chrom = pair.first;
        // iterate over values for each chromosome
        for (unsigned int i = 0; i < pair.second.size(); i++)
        {
            double sum = 0;
            // iterate over all positions i across samples j
            for (unsigned int j = 0; j < all_per_base_coverages.size(); j++)
            {
                // all_per_base_coverages[sample_nr][chromosome][base_pair_position]
                sum += all_per_base_coverages[j][chrom][i]; //sample j, chromosome chrom, position i
            }
            mean_coverage[chrom].push_back(sum / all_per_base_coverages.size()); // (sum over coverage at bp i) / (#nof samples)
        }
    }

}

// returns true if (1 - tolerance) * 10 <= bp_coverage <= (1 + tolerance) * 10 == 8 <= bp_coverage <= 12 for tolerance = 0.2
// bool Averager::in_interval(double current_avg, double bp_coverage)
// {
//     // tolerance of +/- 20 %
//     return bp_coverage >= 0.8 * current_avg && bp_coverage <= 1.2 * current_avg;
// }

// identify ERs with coverage > threshold and length > min_length bp
void Averager::find_ERs(double threshold, int min_length)
{

    int start = 0;
    double current_avg = 0;
    int count = 0;

    //if length > 5 and coverage at bp < 5 -> ER has ended, append
    // if coverage at bp > 5 --> add to current ER
    // if coverage < 5 && length < 5: don't count as ER, reset start and average

    for (auto& pair : mean_coverage)
    {
        for (unsigned int i = 0; i < pair.second.size(); i++)
        {
            double coverage = pair.second[i];
            // coverage is less than threshold
            if (pair.second[i] <= threshold)
            {
                // region at least 5 bp long, append to results
                if ((i - start) > min_length)
                {
                    current_avg /= (i - 1 - start);
                    BedGraphRow expressed_region = BedGraphRow(pair.first, start, i - 1, current_avg);

                    //expressed_region.print();
                    expressed_regions.push_back(expressed_region);

                }
                //region too short to be appended, reset start and current avg
                start = i + 1;
                current_avg = 0;

            }
            // add to current ER
            else if (coverage > threshold)
            {
                current_avg += coverage;
                ++count;
            }
        }
    }
    //std::cout << count << " positions" << std::endl;
}

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
#include <future>



//compute overall mean coverage
void Averager::compute_mean_coverage(std::vector<std::unordered_map<std::string, std::vector<double>>>& all_per_base_coverages, const std::vector<std::string>& chromosome_sequence)
{

    // define workers --> use one thread per chromosome

    std::vector<std::future<std::vector<double>>> workers;
    workers.reserve(chromosome_sequence.size()); //pre-allocate space for workers

    std::cout << "#samples = " << all_per_base_coverages.size() << ", #chromosomes = "<< all_per_base_coverages[0].size() << std::endl;

    // parallel iteration over chromosomes. pair.first = chromosome, pair.second = vector of per base coverages of that chromosome
    for (auto& chrom : chromosome_sequence)
    {
        workers.push_back(std::async(std::launch::async, [&, chrom]{
        std::cout << "COMPUTING MEAN FOR CHROMOSOME " << chrom << std::endl;
        std::vector<double> chrom_mean_vector;
        // iterate over values for each chromosome
        for (unsigned int i = 0; i < all_per_base_coverages[0][chrom].size(); i++)
        {
            double sum = 0.0;
            // iterate over all positions i across samples j
            for (unsigned int j = 0; j < all_per_base_coverages.size(); j++)
            {
                // all_per_base_coverages[sample_nr][chromosome][base_pair_position]
                sum += all_per_base_coverages[j][chrom][i]; //sample j, chromosome chrom, position i
            }
            chrom_mean_vector.push_back(sum / all_per_base_coverages.size()); // (sum over coverage at bp i) / (#nof samples)
        }
            return chrom_mean_vector;
        }));

    }
    //avoid using mutexes, just do single-threaded merge
    // TODO later use >1 thread per chromosome
    for (int k = 0; k < workers.size(); k++)
    {
        std::string chrom = chromosome_sequence[k];
        mean_coverage[chrom] = workers[k].get(); //get result
        assert(all_per_base_coverages[0].at(chrom).size() == mean_coverage[chrom].size());
    }

}


// identify ERs with coverage > threshold and length > min_length bp
void Averager::find_ERs(double threshold, int min_length, const std::vector<std::string>& chromosome_sequence)
{

    int start = 0;
    double current_avg = 0;
    int count = 0;

    //if length > 5 and coverage at bp < 5 -> ER has ended, append
    // if coverage at bp > 5 --> add to current ER
    // if coverage < 5 && length < 5: don't count as ER, reset start and average
    //pair.first = chromosome, pair.second = vector of per base coverages of that chromosome
    for (auto& chrom : chromosome_sequence)
    {
        for (unsigned int i = 0; i < mean_coverage[chrom].size(); i++)
        {
            double coverage =  mean_coverage[chrom][i];
            // coverage is less than threshold
            if ( mean_coverage[chrom][i] <= threshold)
            {
                // region at least 5 bp long, append to results
                if ((i - start) > min_length)
                {
                    current_avg /= (i - 1 - start);
                    BedGraphRow expressed_region = BedGraphRow(chrom, start, i - 1, current_avg);

                    //expressed_region.print();
                    expressed_regions[chrom].push_back(expressed_region);

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

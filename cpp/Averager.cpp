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
void Averager::compute_mean_coverage(std::vector<std::unordered_map<std::string, std::vector<double>>>& all_per_base_coverages)
{

    // define workers --> use one thread per chromosome

    std::unordered_map<std::string, std::future<std::vector<double>>> workers;
    // workers.reserve(all_per_base_coverages.size()); //pre-allocate space for workers

    std::cout << "#samples = " << all_per_base_coverages.size() << ", #chromosomes = "<< all_per_base_coverages[0].size() << std::endl;

    // parallel iteration over chromosomes. pair.first = chromosome, pair.second = vector of per base coverages of that chromosome
    for (auto& item : all_per_base_coverages.at(0))
    {
        std::string chrom = item.first;
        workers[chrom] = std::async(std::launch::async, [&, chrom]{
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
        });

    }
    //avoid using mutexes, just do single-threaded merge
    // TODO later use >1 thread per chromosome

    for (auto& pair : all_per_base_coverages.at(0))
    {
        std::string chrom = pair.first;
        mean_coverage[chrom] = workers[chrom].get(); //get result
        std::cout << "FINISHED MEAN COMPUTATION FOR " << chrom << " with #bp = " << all_per_base_coverages[0].at(chrom).size() << std::endl;
        assert(all_per_base_coverages[0].at(chrom).size() == mean_coverage[chrom].size());
    }

}


// identify ERs with coverage > threshold and length > min_length bp
void Averager::find_ERs(double threshold, int min_length)//, std::vector<std::string>& chromosome_sequence)
{
    // define workers --> use one thread per chromosome
    std::unordered_map<std::string, std::future<std::vector<BedGraphRow>>> workers;
    workers.reserve(mean_coverage.size()); //pre-allocate space for workers

    //if length > 5 and coverage at bp < 5 -> ER has ended, append
    // if coverage at bp > 5 --> add to current ER
    // if coverage < 5 && length < 5: don't count as ER, reset start and average
    //pair.first = chromosome, pair.second = vector of per base coverages of that chromosome
    for (auto& coverage : mean_coverage)
    {
        std::string chrom = coverage.first;

        // PARALLELIZED
        workers[chrom] = std::async(std::launch::async, [&, chrom]
        {
            int start = 0;
            double current_sum = 0;
            int count = 0;
            std::vector<BedGraphRow> chrom_expressed_regions;
            for (unsigned int i = 0; i < mean_coverage[chrom].size(); i++)
            {
                double coverage =  mean_coverage[chrom][i];
                // coverage is less than threshold
                if (mean_coverage[chrom][i] <= threshold)
                {
                    // region at least 5 bp long, append to results
                    if ((i - start) > min_length)
                    {
                        double current_avg = current_sum / (i - start);
                        //++count;
                        assert(count == i - start);
                        std::cout << "count: " << count << ", len = " << i - start << std::endl;
                        BedGraphRow expressed_region = BedGraphRow(chrom, start, i, current_avg);

                        //expressed_region.print();
                        chrom_expressed_regions.push_back(expressed_region);

                    }
                    //region too short to be appended, reset start and current avg
                    start = i + 1;
                    count = 0;
                    current_sum = 0;

                }
                // add to current ER
                else if (coverage > threshold)
                {
                    current_sum += coverage;
                    ++count;
                }

            } // end inner for loop
                return chrom_expressed_regions;
            });
        }

        //avoid using mutexes, just do single-threaded merge
        // TODO later use >1 thread per chromosome
        for (auto& pair : mean_coverage)
        {
            std::string chrom = pair.first;
            expressed_regions[chrom] = workers[chrom].get(); //get result
            //std::cout << "FINISHED MEAN COMPUTATION FOR " << chrom << " with #bp = " << all_per_base_coverages[0].at(chrom).size() << std::endl;
            //assert(all_per_base_coverages[0].at(chrom).size() == mean_coverage[chrom].size());
        }

}
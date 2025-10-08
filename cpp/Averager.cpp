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
// overload input operator for SJRow


// fill vector with coverage value per bp (since different bedgraphs have different binning intervals)
void Averager::compute_per_base_coverage(const BedGraphRow& row, std::vector<double>& per_base_coverage)
{
    // row.end is NOT inclusive
    int position = row.end - row.start; //for just one nt, row.start = 22, row.end = 23 -> position = 1
    assert(position >= 0);
    do
    {
        per_base_coverage.push_back(row.coverage);
        position--;
    }
    while (position > 0);

}
//normalize reads to CPM for better comparability between libraries
//arguments: vector containing per-base coverage of one sample
void Averager::normalize(std::vector<double>& per_base_coverage, const unsigned int library_size)
{
    for (double& e : per_base_coverage)
    {
        e = (e / library_size) * 1e6; //normalize to CPM
    }
}

//read in a .bedgraph file and fill the per_base_coverage and data vectors
//return library size of this sample file
int read_file(const std::string filename,
    std::vector<double>& per_base_coverage,
    std::vector<BedGraphRow>& data)
{

    std::cout << filename << std::endl;
    //read in file from path
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file " << filename << std::endl;
    }
    std::string line;
    int library_size = 0;
    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);
        BedGraphRow row;
        iss >> row.chrom >> row.start >> row.end >> row.coverage;
        // calculate total number of reads that map to this bp interval
        row.total_reads += (row.end - row.start) * (row.coverage); //if start = 22, end = 25, coverage = 3 --> (25 - 22) * 3 = 3 * 3 = 9
        library_size += row.total_reads;
        //row.print();
        compute_per_base_coverage(row, per_base_coverage);
        data.push_back(row);
    }

    //return library size
    return library_size;

}

// read rr file
void read_rr(const std::string filename,
    std::vector<SJRow>& sample_rr)
{

    std::cout << filename << std::endl;
    //read in file from path
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file " << filename << std::endl;
    }
    std::string line;

    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);
        SJRow row = {};
        iss >> row;
        if (!row.left_annotated.empty())
        {
            std::cout << row << std::endl;
        }

        sample_rr.push_back(row);
    }


}


//compute overall average coverage
std::vector<double> Averager::compute_avg_coverage()
{
    std::vector<double> avg_coverage(all_per_base_coverages[0].size());
    //iterate
    std::cout << all_per_base_coverages.size() << "outer dim, inner dim = "<< all_per_base_coverages[1].size() << std::endl;
    //outer loop iterates over each position in each sample
    for (unsigned int i = 0; i < all_per_base_coverages[0].size(); i++)
    {
        double sum = 0;
        // inner loop iterates over each sample to get average of one position across samples
        for (unsigned int j = 0; j < all_per_base_coverages.size(); j++)
        {
            sum += all_per_base_coverages[j][i]; //sample j, position i
        }
        avg_coverage[i] = sum / all_per_base_coverages.size(); // (#nof reads at bp i) / (#nof samples)
    }
    return avg_coverage;
}

// returns true if (1 - tolerance) * 10 <= bp_coverage <= (1 + tolerance) * 10 == 8 <= bp_coverage <= 12 for tolerance = 0.2
bool Averager::in_interval(double current_avg, double bp_coverage)
{
    // tolerance of +/- 20 %
    return bp_coverage >= 0.8 * current_avg && bp_coverage <= 1.2 * current_avg;
}

// identify ERs with coverage > 0.25 and length > 5 bp
std::vector<BedGraphRow> Averager::find_ERs(const std::vector<double>& avg_coverage)
{

    std::vector<BedGraphRow> results;
    int start = 0;
    double current_avg = 0;
    //int end = 1;

    //if length > 5 and coverage at bp < 5 -> ER has ended, append
    // if coverage at bp > 5 --> add to current ER
    // if coverage < 5 && length < 5: don't count as ER, reset start and average

    for (int i = 0; i < avg_coverage.size(); i++)
    {
        // coverage of less than 5
        if (avg_coverage[i] <=0.25 )
        {
            // region at least 5 bp long, append to results
            if ((i - start) > 5)
            {
                current_avg /= (i - 1 - start);
                BedGraphRow der = {"chr19", start, i - 1, current_avg};
                der.length = (i - 1 - start);
                //der.print();
                results.push_back(der);

            }
            //region too short to be appended, reset start and current avg
            start = i + 1;
            current_avg = 0;

        }
        // add to current ER
        else if (avg_coverage[i] > 0.25)
        {
            current_avg += avg_coverage[i];
        }
    }
    return results;
}

int main() {
    std::string bigwig_path = "../data/preprocessing";
    std::string rr_path = "../data/splice_junctions/gtex.junctions.BRAIN.ALL.RR";
    //std::cin >> path; //for later

    // matrix consisting of multiple vectors, where each row in a vector is of type BedGraphRow
    // bedgraphs.size() == 16
    // bedgraphs[0].size() == length of first bedgraph
    std::vector<std::vector<BedGraphRow>> all_bedgraphs;
    std::vector<std::vector<double>> all_per_base_coverages;

    //read in every sample file and create per-base coverage file
    for (const auto & entry : std::filesystem::directory_iterator(bigwig_path))
    {
        std::vector<BedGraphRow> bedgraph;
        std::vector<double> per_base_coverage;
        int library_size = read_file(entry.path().string(), per_base_coverage, bedgraph);

        // normalize read counts
        normalize(per_base_coverage, library_size);

        // add to matrix
        all_bedgraphs.push_back(bedgraph);
        all_per_base_coverages.push_back(per_base_coverage);
    }

    // for each study, group the MM file by sample id and store in unordered map (faster lookup, and order is irrelevant)
    std::unordered_map<std::string, std::vector<BedGraphRow>> sjs_per_sample; //for each sample id, store start, end
    std::vector<SJRow> sample_rr;

    read_rr(rr_path, sample_rr);




    // compute average coverage per read
    std::vector<double> avg_coverage = compute_avg_coverage(all_per_base_coverages);
    //int count = 0;
    // for (int i = 0; i < avg_coverage.size(); ++i)
    // {
    //     if (avg_coverage[i] > 0.25)
    //     {
    //         std::cout << avg_coverage[i] << std::endl;
    //         ++count;
    //     }
    //
    //
    // }
    // std::cout << count << std::endl; // 795778 positions have activity > 0.25

    // find DERs
    std::vector<BedGraphRow> results = find_ERs(avg_coverage);

    std::cout << results.size() << std::endl;
    std::cout << "here" << std::endl;
    // for (auto & result : results)
    // {
    //     result.print();
    // }

    //check dimensions with print statements
    std::cout << all_bedgraphs.size() << std::endl;
    std::cout << all_bedgraphs[0].size()<< std::endl;
    std::cout << all_bedgraphs[1].size()<< std::endl;
    // for (int i = 0; i < 10; ++i)
    // {
    //     for (int j = 0; j < 10; ++j)
    //     {
    //         std::cout << all_per_base_coverages[i][j] << std::endl;
    //     }
    // }

    std::cout << all_per_base_coverages.size() << std::endl;
    std::cout << all_per_base_coverages[0].size() << std::endl;
    std::cout << all_per_base_coverages[1].size() << std::endl;


    return 0;
}
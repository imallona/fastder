//
// Created by martinalavanya on 24.09.25.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>

// custom struct for BedGraph
struct BedGraphRow
{
    std::string chrom;
    int start;
    int end;
    double coverage;
    int total_reads = 0;
    // add optional values for average coverage, DER identifier
};

// fill vector with coverage value per bp (since different bedgraphs have different binning intervals)
void compute_per_base_coverage(const BedGraphRow row, std::vector<double>& per_base_coverage)
{
    // row.end is NOT inclusive
    int position = row.end - row.start; //for just one nt, row.start = 22, row.end = 23 -> position = 1
    do
    {
        per_base_coverage.push_back(row.coverage);
        position--;
    }
    while (position > 0);

}
//normalize reads to CPM for better comparability between libraries
void normalize(std::vector<double>& per_base_coverage, double total_reads)
{
    for (double& e : per_base_coverage)
    {
        e = (e / total_reads) * 1e6; //normalize to CPM
    }
}

//read in a .bedgraph file and fill the per_base_coverage and data vectors
void read_file(const std::string filename,
    std::vector<double>& per_base_coverage,
    std::vector<BedGraphRow>& data)
{

    std::cout << filename << std::endl;
    //read in file from path
    std::ifstream file(filename);
    std::string line;
    // iterate over lines
    while (std::getline(file, line))
    {
        // read in line by line
        std::istringstream iss(line);
        BedGraphRow row;
        iss >> row.chrom >> row.start >> row.end >> row.coverage;
        row.total_reads += (row.end - row.start) * (row.coverage); //if start = 22, end = 25, coverage = 3 --> (25 - 22) * 3 = 3 * 3 = 9
        compute_per_base_coverage(row, per_base_coverage);
        data.push_back(row);

        // normalize read counts
        normalize(per_base_coverage, row.total_reads);
    }

}

//compute overall average coverage
std::vector<double> compute_avg_coverage(const std::vector<std::vector<double>>& all_per_base_coverages)
{
    std::vector<double> avg_coverage(all_per_base_coverages[0].size());
    for (int i = 0; i < all_per_base_coverages.size(); i++)
    {
        double sum = 0;
        for (const double j : all_per_base_coverages[i])
        {
            sum += j;
        }
        avg_coverage[i] = sum / all_per_base_coverages.size(); // (#nof reads at bp i) / (#nof samples)
    }
    return avg_coverage;
}

int main() {
    std::string path = "../data/preprocessing";
    //std::cin >> path; //for later

    // matrix consisting of multiple vectors, where each row in a vector is of type BedGraphRow
    // bedgraphs.size() == 16
    // bedgraphs[0].size() == length of first bedgraph
    std::vector<std::vector<BedGraphRow>> all_bedgraphs;
    std::vector<std::vector<double>> all_per_base_coverages;

    //read in files
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        std::vector<BedGraphRow> bedgraph;
        std::vector<double> per_base_coverage;
        read_file(entry.path().string(), per_base_coverage, bedgraph);
        all_bedgraphs.push_back(bedgraph);
        all_per_base_coverages.push_back(per_base_coverage);
    }
    // compute average coverage per read
    std::vector<double> avg_coverage = compute_avg_coverage(all_per_base_coverages);

    // find DERs

    //check dimensions with print statements
    std::cout << all_bedgraphs.size() << std::endl;
    std::cout << all_bedgraphs[0].size()<< std::endl;
    std::cout << all_bedgraphs[1].size()<< std::endl;
    std::cout << all_per_base_coverages.size() << std::endl;
    std::cout << all_per_base_coverages[0].size() << std::endl;
    std::cout << all_per_base_coverages[1].size() << std::endl;


    return 0;
}
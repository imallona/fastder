//
// Created by martinalavanya on 17.09.25.
//
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

// custom struct for BedGraph
struct BedGraphRowTest
{
    std::string chrom;
    int start;
    int end;
    int coverage;
};

int main2() {
    std::cout << "HELLO" << std::endl;

    //read in file
    std::ifstream file("../gtex.base_sums.BRAIN_GTEX-12ZZX-2826-SM-5BC6K.1.ALL_chr19.bedGraph");
    if (file.is_open())
    {
        std::vector<BedGraphRowTest> data;
        std::vector<int> per_base_coverage;
        // read in the whole file
        std::string line;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            BedGraphRowTest row;
            iss >> row.chrom >> row.start >> row.end >> row.coverage;

            // if multiple nt are binned as one
            int position = row.end - row.start; //for just one nt, row.start = 22, row.end = 23 -> position = 1
            do
            {
                per_base_coverage.push_back(row.coverage);
                position--;
            }
            while (position > 0);
            data.push_back(row);


        }


        std::cout << data.size() << std::endl;
        std::cout << data.back().end << std::endl;
        std::cout << per_base_coverage.size() << std::endl;
        std::cout << data.at(0).chrom << " " << data.at(0).start<< " " << data.at(0).end << " " << data.at(0).coverage;
    }
    for (int i = 0; i < 5; ++i)
    {
        std::vector<int> test;
        test.push_back(i);
        std::cout << test.size() << std::endl;
    }


    return 0;
}
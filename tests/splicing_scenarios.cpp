//
// Created by martinalavanya on 05.11.25.
//

#include "BedGraphRow.h"
#include <gtest/gtest.h>
#include <vector>
#include <unordered_map>
#include <string>

#include "Integrator.h"
#include "SJRow.h"


TEST(SpliceTest, TwoStitchedERsTwoSJs)
{
    // create splice junctions
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false, "CT", "AC", "0", "0"), // id 2
    };

    // create expressed regions map
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100),  // length 500
        BedGraphRow("chr1", 11000, 12500, 101),  // length 1500
        BedGraphRow("chr1", 12861, 12999, 29), // length 138
        BedGraphRow("chr1", 14001, 14540, 30)  // length 539
    };
    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2};
    Integrator integrator = Integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test1.gtf");
    EXPECT_EQ(integrator.stitched_ERs.size(), 2); // two ERs
    EXPECT_EQ(integrator.stitched_ERs.at(0).er_ids.size(), 3); // two ERs, one spliced region
}


TEST(SpliceTest, StitchedERWithThreeERsTwoSJs)
{
    // create splice junctions
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false, "CT", "AC", "0", "0"), // id 2
        SJRow("chr1", 14200, 15000, 800, '-', false, "CT", "AC", "0", "0")  // id 3
    };

    // create expressed regions map
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100),
        BedGraphRow("chr1", 11000, 12500, 101),
        BedGraphRow("chr1", 12861, 12999, 29),
        BedGraphRow("chr1", 14001, 14201, 30),
        BedGraphRow("chr1", 14999, 15300, 30),
    };
    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2, 3};
    Integrator integrator = Integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test2.gtf");
    EXPECT_EQ(integrator.stitched_ERs.size(), 2); // two ERs
    EXPECT_EQ(integrator.stitched_ERs.at(1).er_ids.size(), 5); // three ERs, two spliced region
}

TEST(SpliceTest, StitchedERWithThreeERsTwoSJsAndTailingER)
{
    // create splice junctions
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false, "CT", "AC", "0", "0"), // id 2
        SJRow("chr1", 14200, 15000, 800, '-', false, "CT", "AC", "0", "0")  // id 3
    };

    // create expressed regions map
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100),
        BedGraphRow("chr1", 11000, 12500, 101),
        BedGraphRow("chr1", 12861, 12999, 29),
        BedGraphRow("chr1", 14001, 14201, 30),
        BedGraphRow("chr1", 14999, 15300, 30),
        BedGraphRow("chr1", 15300, 15600, 37),
    };
    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2, 3};
    Integrator integrator = Integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test3.gtf");
    EXPECT_EQ(integrator.stitched_ERs.size(), 3); // three ERs
    EXPECT_EQ(integrator.stitched_ERs.at(2).er_ids.size(), 1);
}

TEST(SpliceTest, StitchedERWithThreeERsTwoSJsAndTwoTailingERs)
{
    // create splice junctions
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false, "CT", "AC", "0", "0"), // id 2
        SJRow("chr1", 14200, 15000, 800, '-', false, "CT", "AC", "0", "0")  // id 3
    };

    // create expressed regions map
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100),
        BedGraphRow("chr1", 11000, 12500, 101),
        BedGraphRow("chr1", 12861, 12999, 29),
        BedGraphRow("chr1", 14001, 14201, 30),
        BedGraphRow("chr1", 14999, 15300, 30),
        BedGraphRow("chr1", 15300, 15600, 37),
        BedGraphRow("chr1", 15600, 15800, 87),
    };
    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2, 3};
    Integrator integrator = Integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test4.gtf");
    EXPECT_EQ(integrator.stitched_ERs.size(), 4); // three ERs
    EXPECT_EQ(integrator.stitched_ERs.at(3).er_ids.size(), 1);
}
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
#include "Parser.h"
#include "Averager.h"

TEST(SpliceTestChromOne, TwoStitchedERsTwoSJs)
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


TEST(SpliceTestChromOne, StitchedERWithThreeERsTwoSJs)
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

TEST(SpliceTestChromOne, StitchedERWithThreeERsTwoSJsAndTailingER)
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

TEST(SpliceTestChromOne, StitchedERWithThreeERsTwoSJsAndTwoTailingERs)
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


TEST(SpliceTestChromOne, NoSJUsed)
{
    // create splice junctions
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 500, 1000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 12000, 13000, 1000, '-', false, "CT", "AC", "0", "0"), // id 2
        SJRow("chr1", 15300, 15400, 100, '-', false, "CT", "AC", "0", "0")  // id 3
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
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test5.gtf");
    EXPECT_EQ(integrator.stitched_ERs.size(), 7); // 7 ERs, no stitching

}


TEST(SpliceTestChromOne, MiddleSJUnusedTailingSJsUsed)
{
    // create splice junctions
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 12000, 13000, 1000, '-', false, "CT", "AC", "0", "0"), // id 2
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
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test6.gtf");

    EXPECT_EQ(integrator.stitched_ERs.size(), 5); // five ERs, 2x2 stitching
    EXPECT_EQ(integrator.stitched_ERs.at(0).er_ids.size(), 3);
    EXPECT_EQ(integrator.stitched_ERs.at(1).er_ids.size(), 1);
    EXPECT_EQ(integrator.stitched_ERs.at(2).er_ids.size(), 3);
    EXPECT_EQ(integrator.stitched_ERs.at(3).er_ids.size(), 1);

}


TEST(SpliceTestChromOne, FirstERNotStitchedRemainingERsStitched)
{
    // create splice junctions
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 14200, 15000, 800, '-', false, "CT", "AC", "0", "0")  // id 3
    };

    // create expressed regions map
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 9000, 9200, 2),
        BedGraphRow("chr1", 10000, 10500, 100),
        BedGraphRow("chr1", 11000, 12500, 101),
        BedGraphRow("chr1", 14001, 14201, 30),
        BedGraphRow("chr1", 14999, 15300, 30),
    };
    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2};
    Integrator integrator = Integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test7.gtf");

    EXPECT_EQ(integrator.stitched_ERs.size(), 3); // five ERs, 2x2 stitching
    EXPECT_EQ(integrator.stitched_ERs.at(0).er_ids.size(), 1);
    EXPECT_EQ(integrator.stitched_ERs.at(1).er_ids.size(), 3);
    EXPECT_EQ(integrator.stitched_ERs.at(2).er_ids.size(), 3);


}

TEST(SpliceTestChromOneAndTwo, BothChrHaveMatchingSJs)
{
    // create splice junctions (ordered within a chromosome but NOT ordered by chromosome!
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr2", 1000, 1500,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr2", 3000, 3800, 800, '-', false, "CT", "AC", "0", "0"),  // id 3
        SJRow("chr1", 10500, 11000,  500, '-', false, "CT", "AC", "0", "0"), // id 1
        SJRow("chr1", 14200, 15000, 800, '-', false, "CT", "AC", "0", "0")  // id 3
    };

    // create expressed regions map
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr2"] = {
        BedGraphRow("chr2", 1501, 2999, 10), // er
        BedGraphRow("chr2", 3798, 4000, 11),
    };
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 9000, 9200, 2), // er
        BedGraphRow("chr1", 10000, 10500, 30), // er
        BedGraphRow("chr1", 11000, 14201, 31),
        BedGraphRow("chr1", 14999, 15300, 32),
    };

    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj;
    mm_chrom_sj["chr2"] = {1, 2};
    mm_chrom_sj["chr1"] = {3, 4};
    Integrator integrator = Integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);
    integrator.write_to_gtf("../../tests/gtfs/splicing_scenarios_test8.gtf");

    int count_chr1 = 0;
    int largest_ser_chr1 = 0;
    int count_chr2 = 0;
    int largest_ser_chr2 = 0;
    for (auto ser : integrator.stitched_ERs)
    {
        if (ser.chrom == "chr1")
        {
            count_chr1++;
            if (ser.er_ids.size() > largest_ser_chr1)
            {
                largest_ser_chr1 = ser.er_ids.size();
            }
        }
        if (ser.chrom == "chr2")
        {
            count_chr2++;
            if (ser.er_ids.size() > largest_ser_chr2)
            {
                largest_ser_chr2 = ser.er_ids.size();
            }
        }
    }
    EXPECT_EQ(integrator.stitched_ERs.size(), 3); // four ERs, 1 on chr2 and 3 on chr1

    EXPECT_EQ(count_chr1, 2);
    EXPECT_EQ(count_chr2, 1);

    EXPECT_EQ(largest_ser_chr1, 5);
    EXPECT_EQ(largest_ser_chr2, 3);
    // EXPECT_EQ(integrator.stitched_ERs.at(1).er_ids.size(), 3);
    // EXPECT_EQ(integrator.stitched_ERs.at(2).er_ids.size(), 3);


}

// TEST(Parser, TestFullWithDummyData)
// {
//     // parse files
//     std::string directory = "../data/test_exon_skipping";
//
//     int position_tolerance = 5;
//     double coverage_tolerance = 0.1;
//
//     // parse files
//     Parser parser = Parser(directory, {});
//     parser.search_directory();
//
//     // get mean coverage vector
//     Averager averager;
//     averager.compute_mean_coverage(parser.all_per_base_coverages);
//
//     // get expressed regions
//     averager.find_ERs(0.25, 5);
//
//     // use splice junctions to stitch together expressed regions
//     Integrator integrator = Integrator(coverage_tolerance, position_tolerance);
//     integrator.stitch_up(averager.expressed_regions, parser.mm_chrom_sj, parser.rr_all_sj);
//
//     // convert to GTF format
//     std::string output_path = "../data/output.gtf";
//     integrator.write_to_gtf(output_path);
// }

TEST(Parser, TestWrongChromosomeOrder)
{
    // parse files
    std::string directory = "../../tests/test_data";

    int position_tolerance = 5;
    int length_threshold = 5;
    double coverage_tolerance = 0.1;
    double coverage_threshold = 0.25;
    std::vector<std::string> chromosomes = {"chr2", "chr1"}; // intentionally wrong order
    // parse files
    Parser parser = Parser(directory, chromosomes);
    parser.search_directory();

    // get mean coverage vector
    Averager averager;
    averager.compute_mean_coverage(parser.all_per_base_coverages);

    // get expressed regions
    averager.find_ERs(coverage_threshold, length_threshold);

    // use splice junctions to stitch together expressed regions
    Integrator integrator = Integrator(coverage_tolerance, position_tolerance);
    integrator.stitch_up(averager.expressed_regions, parser.mm_chrom_sj, parser.rr_all_sj);

    // convert to GTF format
    std::string output_path = "../../tests/gtfs/parser_test2.gtf";
    integrator.write_to_gtf(output_path);
}
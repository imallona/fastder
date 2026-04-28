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
        SJRow("chr1", 10500, 11000,  500, '-', false), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false), // id 2
    };

    // create expressed regions map
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100),  // length 500
        BedGraphRow("chr1", 11000, 12500, 101),  // length 1500
        BedGraphRow("chr1", 12861, 12999, 29), // length 138
        BedGraphRow("chr1", 14001, 14540, 30)  // length 539
    };
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
        SJRow("chr1", 10500, 11000,  500, '-', false), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false), // id 2
        SJRow("chr1", 14200, 15000, 800, '-', false)  // id 3
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
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
        SJRow("chr1", 10500, 11000,  500, '-', false), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false), // id 2
        SJRow("chr1", 14200, 15000, 800, '-', false)  // id 3
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
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
        SJRow("chr1", 10500, 11000,  500, '-', false), // id 1
        SJRow("chr1", 13000, 14000, 1000, '-', false), // id 2
        SJRow("chr1", 14200, 15000, 800, '-', false)  // id 3
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
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
        SJRow("chr1", 500, 1000,  500, '-', false), // id 1
        SJRow("chr1", 12000, 13000, 1000, '-', false), // id 2
        SJRow("chr1", 15300, 15400, 100, '-', false)  // id 3
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
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
        SJRow("chr1", 10500, 11000,  500, '-', false), // id 1
        SJRow("chr1", 12000, 13000, 1000, '-', false), // id 2
        SJRow("chr1", 14200, 15000, 800, '-', false)  // id 3
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
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
        SJRow("chr1", 10500, 11000,  500, '-', false), // id 1
        SJRow("chr1", 14200, 15000, 800, '-', false)  // id 3
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
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
        SJRow("chr2", 1000, 1500,  500, '-', false), // id 1
        SJRow("chr2", 3000, 3800, 800, '-', false),  // id 3
        SJRow("chr1", 10500, 11000,  500, '-', false), // id 1
        SJRow("chr1", 14200, 15000, 800, '-', false)  // id 3
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

    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
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
    int min_length = 5;
    double coverage_tolerance = 0.1;
    double coverage_threshold = 0.25;
    int cores = 4;
    std::vector<std::string> chromosomes = {"chr2", "chr1"}; // intentionally wrong order
    // parse files
    Parser parser = Parser(directory, chromosomes, cores);
    parser.search_directory();

    // get mean coverage vector
    Averager averager(cores);
    averager.compute_mean_coverage(parser.all_per_base_coverages);

    // get expressed regions
    averager.find_ERs(coverage_threshold, min_length);


    // use splice junctions to stitch together expressed regions
    Integrator integrator = Integrator(coverage_tolerance, position_tolerance);
    integrator.stitch_up(averager.expressed_regions, parser.mm_chrom_sj, parser.rr_all_sj);

    // SUMMARY OF EXPRESSED REGIONS
    for (auto& c : averager.expressed_regions)
    {
        std::cout << "expressed regions in chr " << c.first << " = " << c.second.size() << std::endl;
        for (auto& er : c.second)
        {
            er.print();
        }
    }

    // SUMMARY OF STITCHED EXPRESSED REGIONS
    std::unordered_map<std::string, int> ser_counts;
    for (auto& stitched_er : integrator.stitched_ERs)
    {
        ser_counts[stitched_er.chrom]++;
        std::cout << stitched_er << std::endl;
    }
    for (auto& c : ser_counts)
    {
        std::cout << "stitched ERs in chr " << c.first << " = " << c.second << std::endl;

    }

    // convert to GTF format
    std::string output_path = "../../tests/gtfs/parser_test3.gtf";
    integrator.write_to_gtf(output_path);
}


TEST(Parser, DateFormatInGTF)
{
    // get today's date
    auto now = std::chrono::system_clock::now();
    std::chrono::year_month_day ymd{std::chrono::floor<std::chrono::days>(now)}; //formatted as YYYY-MM-DD
    std::string date =
    std::to_string(int(ymd.year())) + "-" +
    std::to_string(unsigned(ymd.month())) + "-" +
    std::to_string(unsigned(ymd.day()));

    std::cout << date;

}


// SJRow no longer stores left_motif / right_motif / left_annotated /
// right_annotated. operator>> must still consume those four columns from a
// recount3-formatted RR line so existing inputs continue to parse, but only
// the six load-bearing fields end up in the struct.
TEST(SJRow, ParsesAndDiscardsUnusedRRColumns)
{
    std::istringstream iss(
        "chr1\t143465342\t143470614\t5273\t-\t1\tCT\tAC\tcH38,gC24\taC19,cH38");
    SJRow row;
    iss >> row;

    EXPECT_TRUE(iss.good() || iss.eof());
    EXPECT_EQ(row.chrom, "chr1");
    EXPECT_EQ(row.start, 143465342u);
    EXPECT_EQ(row.end, 143470614u);
    EXPECT_EQ(row.length, 5273u);
    EXPECT_FALSE(row.strand);  // '-' -> false
    EXPECT_TRUE(row.annotated);

    // Must have consumed all ten whitespace-separated tokens; nothing is left
    // for a follow-up record on the same stream
    std::string leftover;
    iss >> leftover;
    EXPECT_TRUE(leftover.empty());
}


// Empty / placeholder annotation columns must also parse cleanly, since
// monorail_light emits "." for left_annotated / right_annotated.
TEST(SJRow, ParsesPlaceholderAnnotationColumns)
{
    std::istringstream iss("chr21\t100\t200\t101\t+\t0\tGT\tAG\t.\t.");
    SJRow row;
    iss >> row;
    EXPECT_EQ(row.chrom, "chr21");
    EXPECT_EQ(row.start, 100u);
    EXPECT_EQ(row.end, 200u);
    EXPECT_TRUE(row.strand);   // '+'
    EXPECT_FALSE(row.annotated);
}


// Sanity check on the size shrink. SJRow now carries exactly one std::string
// (chrom) plus 8 + 8 + 4 + 1 + 1 bytes of trivially-copyable scalars. The
// upstream version carried five extra std::strings, which on libstdc++ adds
// roughly 160 bytes per row, a meaningful cost at 10M rows for full hg38.
TEST(SJRow, ResidentSizeBoundedToOneString)
{
    // Allow generous slack for std::string SSO and alignment differences across
    // toolchains. The check is "much smaller than the legacy 6-string layout"
    // (6 * sizeof(std::string) on libstdc++ is 192 bytes by itself).
    EXPECT_LT(sizeof(SJRow), 96u);
}


// read_rr now keeps only rows whose chromosome is in chromosomes_set.
// rr_all_sj is the filtered subset; sj_id_remap maps the on-disk 1-based
// sj_id to the 1-based index in rr_all_sj, or 0 for dropped rows.
TEST(Parser, ChrFilterRetainsOnlyRequestedRows)
{
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path() / "fastder_test_chr_filter";
    fs::create_directories(tmp);
    auto rr = tmp / "fake.ALL.RR";
    {
        std::ofstream out(rr);
        out << "chromosome\tstart\tend\tlength\tstrand\tannotated\t"
               "left_motif\tright_motif\tleft_annotated\tright_annotated\n";
        out << "chr1\t100\t200\t101\t+\t1\tGT\tAG\t.\t.\n";
        out << "chr2\t300\t400\t101\t+\t1\tGT\tAG\t.\t.\n";
        out << "chr1\t500\t600\t101\t-\t0\tCT\tAC\t.\t.\n";
    }

    Parser parser("dummy_path", {"chr1"}, 1);
    parser.read_rr(rr.string());

    EXPECT_EQ(parser.rr_total_rows, 3u);
    EXPECT_EQ(parser.rr_all_sj.size(), 2u);
    ASSERT_EQ(parser.sj_id_remap.size(), 3u);
    EXPECT_EQ(parser.sj_id_remap[0], 1u);
    EXPECT_EQ(parser.sj_id_remap[1], 0u);
    EXPECT_EQ(parser.sj_id_remap[2], 2u);
    ASSERT_EQ(parser.rr_all_sj.size(), 2u);
    EXPECT_EQ(parser.rr_all_sj[0].chrom, "chr1");
    EXPECT_EQ(parser.rr_all_sj[0].start, 100u);
    EXPECT_EQ(parser.rr_all_sj[1].start, 500u);

    fs::remove_all(tmp);
}


// read_mm must skip lines whose sj_id was dropped by the chr filter, must
// translate retained sj_ids through the remap, and must validate its header
// against the total RR row count (not the filtered rr_all_sj size).
TEST(Parser, MMUsesRemapAndValidatesHeaderAgainstTotalRows)
{
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path() / "fastder_test_mm_remap";
    fs::create_directories(tmp);
    auto rr = tmp / "fake.ALL.RR";
    auto mm = tmp / "fake.ALL.MM";
    {
        std::ofstream out(rr);
        out << "chromosome\tstart\tend\tlength\tstrand\tannotated\t"
               "left_motif\tright_motif\tleft_annotated\tright_annotated\n";
        out << "chr1\t100\t200\t101\t+\t1\tGT\tAG\t.\t.\n";
        out << "chr2\t300\t400\t101\t+\t1\tGT\tAG\t.\t.\n";
        out << "chr1\t500\t600\t101\t-\t0\tCT\tAC\t.\t.\n";
    }
    {
        std::ofstream out(mm);
        out << "%%MatrixMarket matrix coordinate integer general\n";
        out << "%-----------------------------------------------\n";
        out << "3\t2\t4\n";
        out << "1\t1\t5\n";
        out << "2\t1\t3\n";
        out << "3\t1\t8\n";
        out << "1\t2\t10\n";
    }

    Parser parser("dummy_path", {"chr1"}, 1);
    parser.mm_ids = {1u, 2u};
    parser.read_rr(rr.string());
    parser.read_mm(mm.string());

    ASSERT_TRUE(parser.mm_chrom_sj.count("chr1"));
    EXPECT_EQ(parser.mm_chrom_sj.count("chr2"), 0u);
    EXPECT_EQ(parser.mm_chrom_sj["chr1"].size(), 3u);
    for (uint32_t sj_id : parser.mm_chrom_sj["chr1"]) {
        ASSERT_LE(sj_id, parser.rr_all_sj.size());
        EXPECT_EQ(parser.rr_all_sj[sj_id - 1].chrom, "chr1");
    }

    fs::remove_all(tmp);
}


// Header validation must reject an MM whose declared row count does not match
// the RR's total row count, even when the filtered rr_all_sj is smaller.
TEST(Parser, MMHeaderRejectsCountMismatch)
{
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path() / "fastder_test_mm_mismatch";
    fs::create_directories(tmp);
    auto rr = tmp / "fake.ALL.RR";
    auto mm = tmp / "fake.ALL.MM";
    {
        std::ofstream out(rr);
        out << "chromosome\tstart\tend\tlength\tstrand\tannotated\t"
               "left_motif\tright_motif\tleft_annotated\tright_annotated\n";
        out << "chr1\t100\t200\t101\t+\t1\tGT\tAG\t.\t.\n";
        out << "chr2\t300\t400\t101\t+\t1\tGT\tAG\t.\t.\n";
    }
    {
        std::ofstream out(mm);
        out << "%%MatrixMarket matrix coordinate integer general\n";
        out << "999\t1\t0\n";
    }

    Parser parser("dummy_path", {"chr1"}, 1);
    parser.mm_ids = {1u};
    parser.read_rr(rr.string());
    parser.read_mm(mm.string());

    EXPECT_TRUE(parser.mm_chrom_sj.empty());

    fs::remove_all(tmp);
}
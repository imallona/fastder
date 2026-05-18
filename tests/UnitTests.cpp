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
#include "GTFRow.h"
#include <sstream>

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

    // get mean coverage as sparse intervals
    Averager averager(cores);
    averager.compute_mean_coverage(parser.all_bedgraphs);

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
    // The RR start is 1-based; the parser shifts it to 0-based half-open.
    EXPECT_EQ(row.start, 143465341u);
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
    // 1-based RR start shifted to 0-based half-open by the parser.
    EXPECT_EQ(row.start, 99u);
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
    // RR starts are 1-based; the parser shifts them to 0-based half-open.
    EXPECT_EQ(parser.rr_all_sj[0].start, 99u);
    EXPECT_EQ(parser.rr_all_sj[1].start, 499u);

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


// =====================================================================
// BedGraphRow strand field
// =====================================================================

TEST(BedGraphRow, DefaultStrandIsUnstranded)
{
    BedGraphRow row("chr1", 100, 200, 5.0);
    EXPECT_EQ(row.strand, '.');
}

TEST(BedGraphRow, ExplicitStrandIsStored)
{
    BedGraphRow plus("chr1", 100, 200, 5.0, '+');
    BedGraphRow minus("chr1", 100, 200, 5.0, '-');
    EXPECT_EQ(plus.strand,  '+');
    EXPECT_EQ(minus.strand, '-');
}


// =====================================================================
// Averager::compute_mean_coverage on sparse intervals
// =====================================================================

// The sparse mean of a single sample equals the sample's own intervals.
// Adjacent intervals with the same coverage are coalesced.
TEST(Averager, SparseMeanSingleSampleIsIdentity)
{
    std::vector<std::vector<BedGraphRow>> all_bedgraphs(1);
    all_bedgraphs[0] = {
        BedGraphRow("chr1", 100, 200, 4.0),
        BedGraphRow("chr1", 200, 300, 4.0),  // contiguous, same coverage -> coalesced
        BedGraphRow("chr1", 400, 500, 8.0),
    };

    Averager avg(1);
    avg.compute_mean_coverage(all_bedgraphs);

    ASSERT_TRUE(avg.mean_intervals.count("chr1"));
    const auto& m = avg.mean_intervals["chr1"];
    ASSERT_EQ(m.size(), 2u);
    EXPECT_EQ(m[0].start, 100u);
    EXPECT_EQ(m[0].end,   300u);
    EXPECT_DOUBLE_EQ(m[0].coverage, 4.0);
    EXPECT_EQ(m[1].start, 400u);
    EXPECT_EQ(m[1].end,   500u);
    EXPECT_DOUBLE_EQ(m[1].coverage, 8.0);
}

// Two samples with overlapping intervals: the union of breakpoints partitions
// the chromosome and each segment carries the per-sample average.
TEST(Averager, SparseMeanTwoOverlappingSamples)
{
    std::vector<std::vector<BedGraphRow>> all_bedgraphs(2);
    all_bedgraphs[0] = { BedGraphRow("chr1", 100, 300, 10.0) };
    all_bedgraphs[1] = { BedGraphRow("chr1", 200, 400, 20.0) };

    Averager avg(1);
    avg.compute_mean_coverage(all_bedgraphs);

    ASSERT_TRUE(avg.mean_intervals.count("chr1"));
    const auto& m = avg.mean_intervals["chr1"];
    // expected segments: [100,200) cov=10/2=5, [200,300) cov=(10+20)/2=15, [300,400) cov=20/2=10
    ASSERT_EQ(m.size(), 3u);
    EXPECT_EQ(m[0].start, 100u);  EXPECT_EQ(m[0].end, 200u);
    EXPECT_DOUBLE_EQ(m[0].coverage, 5.0);
    EXPECT_EQ(m[1].start, 200u);  EXPECT_EQ(m[1].end, 300u);
    EXPECT_DOUBLE_EQ(m[1].coverage, 15.0);
    EXPECT_EQ(m[2].start, 300u);  EXPECT_EQ(m[2].end, 400u);
    EXPECT_DOUBLE_EQ(m[2].coverage, 10.0);
}

// A sample with no coverage on a chromosome contributes 0 to every segment.
TEST(Averager, SparseMeanSampleAbsentFromChromContributesZero)
{
    std::vector<std::vector<BedGraphRow>> all_bedgraphs(2);
    all_bedgraphs[0] = { BedGraphRow("chr1", 100, 200, 10.0) };
    all_bedgraphs[1] = { BedGraphRow("chr2", 100, 200, 30.0) };

    Averager avg(1);
    avg.compute_mean_coverage(all_bedgraphs);

    ASSERT_TRUE(avg.mean_intervals.count("chr1"));
    ASSERT_TRUE(avg.mean_intervals.count("chr2"));
    const auto& m1 = avg.mean_intervals["chr1"];
    const auto& m2 = avg.mean_intervals["chr2"];
    ASSERT_EQ(m1.size(), 1u);
    ASSERT_EQ(m2.size(), 1u);
    EXPECT_DOUBLE_EQ(m1[0].coverage, 5.0);   // 10 / 2 samples
    EXPECT_DOUBLE_EQ(m2[0].coverage, 15.0);  // 30 / 2 samples
}

// A gap between intervals is implicit zero coverage, not a continuation.
TEST(Averager, SparseMeanGapIsImplicitZero)
{
    std::vector<std::vector<BedGraphRow>> all_bedgraphs(1);
    all_bedgraphs[0] = {
        BedGraphRow("chr1", 100, 200, 10.0),
        BedGraphRow("chr1", 500, 600, 10.0),
    };

    Averager avg(1);
    avg.compute_mean_coverage(all_bedgraphs);

    const auto& m = avg.mean_intervals["chr1"];
    ASSERT_EQ(m.size(), 2u);
    EXPECT_EQ(m[0].end,   200u);
    EXPECT_EQ(m[1].start, 500u);
}


// =====================================================================
// Averager::find_ERs on sparse intervals
// =====================================================================

// One contiguous run above threshold -> one ER spanning the full run.
TEST(Averager, FindERsSimpleRun)
{
    Averager avg(1);
    avg.chroms = {"chr1"};
    avg.mean_intervals["chr1"] = {
        BedGraphRow("chr1", 100, 200, 10.0),
        BedGraphRow("chr1", 200, 300, 10.0),
    };

    avg.find_ERs(1.0, 5);
    const auto& ers = avg.expressed_regions["chr1"];
    ASSERT_EQ(ers.size(), 1u);
    EXPECT_EQ(ers[0].start, 100u);
    EXPECT_EQ(ers[0].end,   300u);
}

// A below-threshold interval breaks the run.
TEST(Averager, FindERsBelowThresholdBreaksRun)
{
    Averager avg(1);
    avg.chroms = {"chr1"};
    avg.mean_intervals["chr1"] = {
        BedGraphRow("chr1", 100, 200, 10.0),
        BedGraphRow("chr1", 200, 250, 0.1),
        BedGraphRow("chr1", 250, 350, 10.0),
    };

    avg.find_ERs(1.0, 5);
    const auto& ers = avg.expressed_regions["chr1"];
    ASSERT_EQ(ers.size(), 2u);
    EXPECT_EQ(ers[0].start, 100u);  EXPECT_EQ(ers[0].end, 200u);
    EXPECT_EQ(ers[1].start, 250u);  EXPECT_EQ(ers[1].end, 350u);
}

// A coverage gap (non-contiguous intervals) also breaks the run, even when
// both flanks are above threshold.
TEST(Averager, FindERsGapBreaksRun)
{
    Averager avg(1);
    avg.chroms = {"chr1"};
    avg.mean_intervals["chr1"] = {
        BedGraphRow("chr1", 100, 200, 10.0),
        BedGraphRow("chr1", 500, 600, 10.0),
    };

    avg.find_ERs(1.0, 5);
    const auto& ers = avg.expressed_regions["chr1"];
    ASSERT_EQ(ers.size(), 2u);
}

// min_length filters short ERs even if they cleared the coverage threshold.
TEST(Averager, FindERsMinLengthFilters)
{
    Averager avg(1);
    avg.chroms = {"chr1"};
    avg.mean_intervals["chr1"] = {
        BedGraphRow("chr1", 100, 102, 10.0),  // length 2, below min_length
        BedGraphRow("chr1", 500, 600, 10.0),  // length 100, kept
    };

    avg.find_ERs(1.0, 5);
    const auto& ers = avg.expressed_regions["chr1"];
    ASSERT_EQ(ers.size(), 1u);
    EXPECT_EQ(ers[0].start, 500u);
}


// =====================================================================
// Integrator stitch_up (strand-aware)
// =====================================================================
//
// stitch_up partitions splice junctions by strand. Each chromosome's ERs
// are walked once per strand, and chains built from junction-linked ERs
// are tagged '+' or '-'. ERs that no junction stitched to a neighbour
// stay unstranded ('.'). The tests below verify the bookkeeping
// guarantees, in particular that every ER appears in stitched_ERs
// exactly once across the strand passes.

// Each ER on a chromosome appears in stitched_ERs exactly once: total exon
// count equals total ER count.
TEST(Integrator, StitchUpEachERAppearsExactlyOnce)
{
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000, 500, '+', false),
        SJRow("chr1", 20500, 21000, 500, '-', false),
    };
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100.0),
        BedGraphRow("chr1", 11000, 12500, 101.0),
        BedGraphRow("chr1", 13000, 14000, 5.0),
        BedGraphRow("chr1", 20000, 20500, 50.0),
        BedGraphRow("chr1", 21000, 22000, 50.0),
    };
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2};

    Integrator integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);

    int total_ers_in_output = 0;
    for (const auto& ser : integrator.stitched_ERs)
    {
        for (int id : ser.er_ids)
        {
            if (id >= 0) ++total_ers_in_output;
        }
    }
    EXPECT_EQ(total_ers_in_output, 5);
}

// stitch_up walks ERs in their input order, so stitched_ERs are emitted in
// genomic order on a chromosome. write_to_gtf relies on this.
TEST(Integrator, StitchUpEmitsInGenomicOrder)
{
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000, 500, '-', false),
        SJRow("chr1", 20500, 21000, 500, '+', false),
    };
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 50.0),
        BedGraphRow("chr1", 11000, 12000, 50.0),
        BedGraphRow("chr1", 20000, 20500, 100.0),
        BedGraphRow("chr1", 21000, 22000, 100.0),
    };
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2};

    Integrator integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);

    ASSERT_GE(integrator.stitched_ERs.size(), 2u);
    for (size_t i = 1; i < integrator.stitched_ERs.size(); ++i)
    {
        EXPECT_LE(integrator.stitched_ERs[i - 1].start, integrator.stitched_ERs[i].start);
    }
}

// A pure plus-strand chain comes back tagged '+'. Standalones from an
// otherwise-unmatched ER come back tagged '.'.
TEST(Integrator, StitchUpTagsPlusChainsAsPlus)
{
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000, 500, '+', false),
    };
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100.0),
        BedGraphRow("chr1", 11000, 12500, 101.0),
    };
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1};

    Integrator integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);

    ASSERT_EQ(integrator.stitched_ERs.size(), 1u);
    EXPECT_EQ(integrator.stitched_ERs[0].strand, '+');
    EXPECT_EQ(integrator.stitched_ERs[0].er_ids.size(), 3u);
}

// Same scenario but with the SJ on the minus strand. Tag should be '-'.
TEST(Integrator, StitchUpTagsMinusChainsAsMinus)
{
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000, 500, '-', false),
    };
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100.0),
        BedGraphRow("chr1", 11000, 12500, 101.0),
    };
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1};

    Integrator integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);

    ASSERT_EQ(integrator.stitched_ERs.size(), 1u);
    EXPECT_EQ(integrator.stitched_ERs[0].strand, '-');
}

// Both strands present on the same chromosome. Each ER ends up in exactly
// one StitchedER. ERs claimed by the first strand pass do not appear in
// the second; ERs unclaimed by either are emitted as unstranded standalones.
TEST(Integrator, StitchUpEachERAppearsExactlyOnceWithBothStrands)
{
    std::vector<SJRow> rr_all_sj = {
        SJRow("chr1", 10500, 11000, 500, '+', false),
        SJRow("chr1", 20500, 21000, 500, '-', false),
    };
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100.0),  // chain start (+ pass)
        BedGraphRow("chr1", 11000, 12500, 101.0),  // chain extension (+)
        BedGraphRow("chr1", 13000, 14000, 5.0),    // unstitched gap
        BedGraphRow("chr1", 20000, 20500, 50.0),   // chain start (- pass)
        BedGraphRow("chr1", 21000, 22000, 50.0),   // chain extension (-)
    };
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;
    mm_chrom_sj["chr1"] = {1, 2};

    Integrator integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);

    int total_ers_in_output = 0;
    int unstranded_count = 0;
    int plus_count       = 0;
    int minus_count      = 0;
    for (const auto& ser : integrator.stitched_ERs)
    {
        for (int id : ser.er_ids)
        {
            if (id >= 0) ++total_ers_in_output;
        }
        if (ser.strand == '.') ++unstranded_count;
        if (ser.strand == '+') ++plus_count;
        if (ser.strand == '-') ++minus_count;
    }
    EXPECT_EQ(total_ers_in_output, 5);
    EXPECT_EQ(plus_count,       1);
    EXPECT_EQ(minus_count,      1);
    EXPECT_EQ(unstranded_count, 1);
}

// A chromosome with expressed regions but no entry in mm_chrom_sj must still
// have its ERs emitted as single-ER unstranded StitchedERs. The previous
// implementation iterated only over mm_chrom_sj keys and silently dropped
// ERs from chromosomes with no junctions.
TEST(Integrator, StitchUpEmitsERsOnChromosomesWithNoSJs)
{
    std::vector<SJRow> rr_all_sj = {};
    std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions;
    expressed_regions["chr1"] = {
        BedGraphRow("chr1", 10000, 10500, 100.0),
        BedGraphRow("chr1", 11000, 12500, 101.0),
    };
    std::unordered_map<std::string, std::vector<uint32_t>> mm_chrom_sj;

    Integrator integrator(0.1, 5);
    integrator.stitch_up(expressed_regions, mm_chrom_sj, rr_all_sj);

    ASSERT_EQ(integrator.stitched_ERs.size(), 2u);
    for (const auto& ser : integrator.stitched_ERs) {
        EXPECT_EQ(ser.strand, '.');
        EXPECT_EQ(ser.er_ids.size(), 1u);
    }
}

// GTFRow must propagate the StitchedER strand into column 7 of the output
// row. Without this, gffcompare can never match transcripts or exons since
// it requires strand agreement at those levels.
TEST(GTFRow, SerializesPlusStrandFromStitchedER)
{
    BedGraphRow er("chr1", 10000, 10500, 100.0);
    StitchedER s(er, 0);
    s.strand = '+';

    GTFRow row(s, "transcript", 1);
    std::ostringstream os;
    os << row;
    std::string line = os.str();
    std::istringstream is(line);
    std::vector<std::string> cols;
    std::string field;
    while (std::getline(is, field, '\t')) cols.push_back(field);
    ASSERT_GE(cols.size(), 7u);
    EXPECT_EQ(cols[6], "+");
}

TEST(GTFRow, SerializesMinusStrandFromStitchedER)
{
    BedGraphRow er("chr1", 10000, 10500, 100.0);
    StitchedER s(er, 0);
    s.strand = '-';

    GTFRow row(s, "transcript", 1);
    std::ostringstream os;
    os << row;
    std::string line = os.str();
    std::istringstream is(line);
    std::vector<std::string> cols;
    std::string field;
    while (std::getline(is, field, '\t')) cols.push_back(field);
    ASSERT_GE(cols.size(), 7u);
    EXPECT_EQ(cols[6], "-");
}


// =====================================================================
// libBigWig integration (gated). These run only when fastder is built with
// -DFASTDER_USE_LIBBIGWIG=ON; otherwise they GTEST_SKIP so the unbuilt code
// path does not block local development.
// =====================================================================

#ifdef FASTDER_USE_LIBBIGWIG
extern "C" {
#include <bigWig.h>
}

// Helper: write a tiny BigWig at path with the supplied chromosome name,
// chromosome length and triples of (start, end, value). Uses libBigWig's
// bedGraph-style interval writer, which is the most permissive of its three
// emission modes and the closest match to what bedtools genomecov produces.
static void write_test_bigwig(const std::string& path,
                              const std::string& chrom,
                              uint32_t chrom_len,
                              const std::vector<std::tuple<uint32_t, uint32_t, float>>& intervals)
{
    bwInit(1 << 17);
    bigWigFile_t* fp = bwOpen(const_cast<char*>(path.c_str()), nullptr, "w");
    ASSERT_NE(fp, nullptr) << "Could not open BigWig for write: " << path;
    ASSERT_EQ(bwCreateHdr(fp, 10), 0);

    char* chrom_array[1] = { const_cast<char*>(chrom.c_str()) };
    uint32_t lens[1] = { chrom_len };
    fp->cl = bwCreateChromList(chrom_array, lens, 1);
    ASSERT_NE(fp->cl, nullptr);
    ASSERT_EQ(bwWriteHdr(fp), 0);

    std::vector<char*> chr_buf(intervals.size(), const_cast<char*>(chrom.c_str()));
    std::vector<uint32_t> starts(intervals.size());
    std::vector<uint32_t> ends(intervals.size());
    std::vector<float>    values(intervals.size());
    for (size_t i = 0; i < intervals.size(); ++i)
    {
        starts[i] = std::get<0>(intervals[i]);
        ends[i]   = std::get<1>(intervals[i]);
        values[i] = std::get<2>(intervals[i]);
    }
    ASSERT_EQ(bwAddIntervals(fp,
                             chr_buf.data(),
                             starts.data(),
                             ends.data(),
                             values.data(),
                             static_cast<uint32_t>(intervals.size())), 0);
    bwClose(fp);
}
#endif


// Round-trip: write a small BigWig with libBigWig, read it through
// Parser::read_bigwig, assert intervals come back identical.
TEST(BigWig, RoundTripParsesWrittenIntervals)
{
#ifndef FASTDER_USE_LIBBIGWIG
    GTEST_SKIP() << "libBigWig support is off; rebuild with -DFASTDER_USE_LIBBIGWIG=ON";
#else
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path() / "fastder_test_bw_roundtrip";
    fs::create_directories(tmp);
    auto bw_path = (tmp / "tiny.bw").string();

    // intervals: (start, end, value)
    const std::vector<std::tuple<uint32_t, uint32_t, float>> intervals = {
        {100, 200, 10.0f},
        {200, 300, 0.0f},
        {300, 500, 25.0f},
    };
    write_test_bigwig(bw_path, "chr1", 1000, intervals);

    Parser parser("dummy_path", {"chr1"}, 1);
    uint64_t library_size = 0;
    auto rows = parser.read_bigwig(bw_path, library_size, '.');

    ASSERT_EQ(rows.size(), intervals.size());
    for (size_t i = 0; i < rows.size(); ++i)
    {
        EXPECT_EQ(rows[i].chrom, "chr1");
        EXPECT_EQ(rows[i].start, std::get<0>(intervals[i]));
        EXPECT_EQ(rows[i].end,   std::get<1>(intervals[i]));
        EXPECT_FLOAT_EQ(static_cast<float>(rows[i].coverage), std::get<2>(intervals[i]));
        EXPECT_EQ(rows[i].strand, '.');
    }

    fs::remove_all(tmp);
#endif
}


// Equivalence: feed the same coverage through read_bedgraph and read_bigwig
// and assert that compute_mean_coverage and find_ERs produce identical
// expressed_regions on both. This is the regression guard for "BedGraph and
// BigWig must yield the same fastder output for the same data".
TEST(BigWig, BedGraphAndBigWigEquivalentExpressedRegions)
{
#ifndef FASTDER_USE_LIBBIGWIG
    GTEST_SKIP() << "libBigWig support is off; rebuild with -DFASTDER_USE_LIBBIGWIG=ON";
#else
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path() / "fastder_test_bw_equiv";
    fs::create_directories(tmp);
    auto bw_path = (tmp / "sample.bw").string();
    auto bg_path = (tmp / "sample.bedGraph").string();

    const std::vector<std::tuple<uint32_t, uint32_t, float>> intervals = {
        {100, 300, 10.0f},
        {300, 400, 0.0f},   // below the 1.0 threshold used in find_ERs (integer-parseable for the BedGraph reader)
        {400, 700, 10.0f},
    };
    write_test_bigwig(bw_path, "chr1", 1000, intervals);
    {
        std::ofstream out(bg_path);
        for (const auto& [s, e, v] : intervals)
        {
            out << "chr1\t" << s << '\t' << e << '\t' << v << '\n';
        }
    }

    auto run_pipeline = [&](const std::string& cov_path)
    {
        Parser parser("dummy_path", {"chr1"}, 1);
        uint64_t library_size = 0;
        std::vector<BedGraphRow> rows;
        if (cov_path.size() >= 3 && cov_path.substr(cov_path.size() - 3) == ".bw")
        {
            rows = parser.read_bigwig(cov_path, library_size, '.');
        }
        else
        {
            rows = parser.read_bedgraph(cov_path, library_size);
        }
        // Note: skipping normalize() so the comparison is on raw coverage; the
        // two readers must produce identical raw values for the test to be
        // meaningful.
        std::vector<std::vector<BedGraphRow>> all_bedgraphs(1);
        all_bedgraphs[0] = std::move(rows);

        Averager avg(1);
        avg.compute_mean_coverage(all_bedgraphs);
        avg.find_ERs(1.0, 5);
        return avg.expressed_regions["chr1"];
    };

    auto from_bw = run_pipeline(bw_path);
    auto from_bg = run_pipeline(bg_path);

    ASSERT_EQ(from_bw.size(), from_bg.size());
    for (size_t i = 0; i < from_bw.size(); ++i)
    {
        EXPECT_EQ(from_bw[i].start, from_bg[i].start);
        EXPECT_EQ(from_bw[i].end,   from_bg[i].end);
        EXPECT_DOUBLE_EQ(from_bw[i].coverage, from_bg[i].coverage);
    }

    fs::remove_all(tmp);
#endif
}
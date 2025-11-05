// //
// // Created by martinalavanya on 05.11.25.
// //
//
// #include <gtest/gtest.h>
// #include <unordered_map>
// #include <map>
// #include <vector>
// #include <string>
// #include <cstdint>
//
// TEST(StitchUp, AlternativeFirstExon) {
//     // chr1: E0a -> E1, E0b -> E1 (two alternative first exons converging)
//     // ERs in order: E0a(50-80), E0b(60-90), E1(300-400)
//     std::unordered_map<std::string, std::vector<BedGraphRow>> expressed{
//         {"chr1", { ER("chr1", 50, 80), ER("chr1", 60, 90), ER("chr1", 300, 400) }}
//     };
//     std::vector<SJRow> rr_all {
//         SJ("chr1", 80, 300), // E0a -> E1
//         SJ("chr1", 90, 300), // E0b -> E1
//     };
//     std::map<std::string, std::vector<uint64_t>> mm {
//         {"chr1", {0, 1}}
//     };
//
//     Integrator I;
//     I.stitch_up(expressed, mm, rr_all);
//
//     // Expect two stitched paths or at least that E1 can be stitched from either start.
//     // Given the algorithm is greedy, we validate that at least one stitched ER chain includes E1 after a first exon.
//     // Check there exists a stitched_ER whose er_ids reflect either [0,2] or [1,2].
//     bool has_0_2 = false, has_1_2 = false;
//     for (auto& se : I.stitched_ERs) {
//         if (se.er_ids.size() == 2 && se.er_ids[0] == 0 && se.er_ids[1] == 2) has_0_2 = true;
//         if (se.er_ids.size() == 2 && se.er_ids[0] == 1 && se.er_ids[1] == 2) has_1_2 = true;
//     }
//     EXPECT_TRUE(has_0_2 || has_1_2);
// }
//
// TEST(StitchUp, AlternativeLastExon) {
//     // chr1: E1 -> E2a and E1 -> E2b (two alternative last exons)
//     // ERs: E1(100-200), E2a(300-380), E2b(300-420)
//     std::unordered_map<std::string, std::vector<BedGraphRow>> expressed{
//         {"chr1", { ER("chr1", 100, 200), ER("chr1", 300, 380), ER("chr1", 300, 420) }}
//     };
//     std::vector<SJRow> rr_all {
//         SJ("chr1", 200, 300), // E1 -> E2a
//         SJ("chr1", 200, 300), // E1 -> E2b (same acceptor start, different exon length)
//     };
//     std::map<std::string, std::vector<uint64_t>> mm {
//         {"chr1", {0, 1}}
//     };
//
//     Integrator I;
//     I.stitch_up(expressed, mm, rr_all);
//
//     // Expect two stitched ERs: [0,1] and [0,2] in some order.
//     bool has_0_1 = false, has_0_2 = false;
//     for (auto& se : I.stitched_ERs) {
//         if (se.er_ids.size() == 2 && se.er_ids[0] == 0 && se.er_ids[1] == 1) has_0_1 = true;
//         if (se.er_ids.size() == 2 && se.er_ids[0] == 0 && se.er_ids[1] == 2) has_0_2 = true;
//     }
//     EXPECT_TRUE(has_0_1);
//     EXPECT_TRUE(has_0_2);
// }
//
// TEST(StitchUp, ExonSkipping_Cassette) {
//     // chr1: E1 -> E2 -> E3, plus a skipping junction E1 -> E3
//     // ERs: E1(100-200), E2(300-350), E3(500-600)
//     std::unordered_map<std::string, std::vector<BedGraphRow>> expressed{
//         {"chr1", { ER("chr1", 100, 200), ER("chr1", 300, 350), ER("chr1", 500, 600) }}
//     };
//     std::vector<SJRow> rr_all {
//         SJ("chr1", 200, 300), // E1 -> E2
//         SJ("chr1", 350, 500), // E2 -> E3
//         SJ("chr1", 200, 500), // E1 -> E3 (skip E2)
//     };
//     // Order SJs so the greedy walk will see both the canonical path and the skipping path
//     std::map<std::string, std::vector<uint64_t>> mm {
//         {"chr1", {0, 1, 2}}
//     };
//
//     Integrator I;
//     I.stitch_up(expressed, mm, rr_all);
//
//     // Expect at least one stitched chain [0,1,2] (E1->E2->E3) OR [0,2] (E1->E3).
//     bool has_0_1_2 = false, has_0_2 = false;
//     for (auto& se : I.stitched_ERs) {
//         if (se.er_ids.size() == 3 && se.er_ids[0] == 0 && se.er_ids[1] == 1 && se.er_ids[2] == 2) has_0_1_2 = true;
//         if (se.er_ids.size() == 2 && se.er_ids[0] == 0 && se.er_ids[1] == 2) has_0_2 = true;
//     }
//     EXPECT_TRUE(has_0_1_2 || has_0_2);
// }
//
// TEST(StitchUp, MultipleExonSkipping) {
//     // chr1: E1 -> E2 -> E3 -> E4, plus skipping junction E1 -> E4 (skip E2,E3)
//     // ERs: E1(100-200), E2(300-350), E3(450-480), E4(600-700)
//     std::unordered_map<std::string, std::vector<BedGraphRow>> expressed{
//         {"chr1", {
//             ER("chr1", 100, 200),
//             ER("chr1", 300, 350),
//             ER("chr1", 450, 480),
//             ER("chr1", 600, 700)
//         }}
//     };
//     std::vector<SJRow> rr_all {
//         SJ("chr1", 200, 300), // E1 -> E2
//         SJ("chr1", 350, 450), // E2 -> E3
//         SJ("chr1", 480, 600), // E3 -> E4
//         SJ("chr1", 200, 600)  // E1 -> E4 (skip E2,E3)
//     };
//     std::map<std::string, std::vector<uint64_t>> mm {
//         {"chr1", {0, 1, 2, 3}}
//     };
//
//     Integrator I;
//     I.stitch_up(expressed, mm, rr_all);
//
//     // Expect at least one stitched chain [0,1,2,3] (full) OR [0,3] (skip two exons).
//     bool has_full = false, has_skip = false;
//     for (auto& se : I.stitched_ERs) {
//         if (se.er_ids.size() == 4
//             && se.er_ids[0] == 0 && se.er_ids[1] == 1 && se.er_ids[2] == 2 && se.er_ids[3] == 3) has_full = true;
//         if (se.er_ids.size() == 2
//             && se.er_ids[0] == 0 && se.er_ids[1] == 3) has_skip = true;
//     }
//     EXPECT_TRUE(has_full || has_skip);
// }
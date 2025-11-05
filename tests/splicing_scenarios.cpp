//
// Created by martinalavanya on 05.11.25.
//

#include "BedGraphRow.h"
#include <gtest/gtest.h>


TEST(BedGraphRowTest, ConstructorInitializesValues)
{
    BedGraphRow row("chr1", 100, 200, 2.5);
    EXPECT_EQ(row.chrom, "chr1");
    EXPECT_EQ(row.start, 100);
    EXPECT_EQ(row.end, 200);
    EXPECT_EQ(row.length, 100);
}
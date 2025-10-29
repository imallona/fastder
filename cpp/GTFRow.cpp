//
// Created by martinalavanya on 28.10.25.
//

#include "GTFRow.h"

#include <iostream>
#include <regex>


GTFRow::GTFRow(const StitchedER& region)
{
    // stitchedER has only one exon
    if (region.er_ids.size() == 1)
    {
        feature = "exon";
    }
    else if (region.er_ids.size() > 1)
    {
        feature = "gene";
    }
    else
    {
        feature = "UNKNOWN";
        std::cout << "ERROR: Stitched ER contains " << region.er_ids.size() << " exons!";
    }
    seqname = std::regex_replace(region.chrom, std::regex("^chr"), ""); //strip away chr from chr1, chrX etc.
    score = region.across_er_coverage;
    start = region.start;
    end = region.end;
    attribute = "nof_expressed_regions="; //hid=trf; hstart=1; hend=21
    attribute += std::to_string(region.er_ids.size());
    attribute += "; length=";
    attribute += std::to_string(region.total_length);

}

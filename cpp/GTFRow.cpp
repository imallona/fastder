//
// Created by martinalavanya on 28.10.25.
//

#include "GTFRow.h"

#include <iostream>
#include <regex>


GTFRow::GTFRow(const StitchedER& region, std::string ftr, unsigned int id)
{

    seqname = std::regex_replace(region.chrom, std::regex("^chr"), ""); //strip away chr from chr1, chrX etc.
    feature = ftr;
    score = region.across_er_coverage;
    start = region.start;
    end = region.end;
    // update column 9 based on feature
    change_feature(ftr, id, 0);


}

void GTFRow::change_feature(std::string ftr, unsigned int id, unsigned int exon_nr)
{
    feature = ftr;
    // always include gene id
    attribute = "gene_id \"gene";
    attribute += std::to_string(id);

    if (ftr == "gene")
    {
        attribute += "\"; gene_name \"faster_gene";
        attribute += std::to_string(id);
    }
    else if (ftr == "transcript" || ftr == "exon")
    {
        attribute += "\"; transcript_id \"tx";
        attribute += std::to_string(id);
    }
    else
    {
        std::cout << "ERROR: UNKNOWN FEATURE";
    }

    // add exon number only if it's an exon
    if (ftr == "exon" && exon_nr > 0)
    {
        attribute += "\"; exon_number \"";
        attribute += std::to_string(exon_nr);
    }


    attribute += "\";";
    // attribute += "\"; nof_expressed_regions="; //hid=trf; hstart=1; hend=21
    // attribute += std::to_string(region.er_ids.size());
    // attribute += "; length=";
    // attribute += std::to_string(region.total_length);

}
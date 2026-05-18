//
// Created by martinalavanya on 20.10.25.
//

#ifndef MLS_STITCHEDER_H
#define MLS_STITCHEDER_H

#include <vector>
#include <BedGraphRow.h>
#include <cstdint>

class StitchedER
{
public:
    // MEMBER FUNCTIONS
    StitchedER() = default;

    StitchedER(const BedGraphRow& expressed_region, int er_id);
    void append(int er_id, unsigned int length, double coverage,
                uint64_t lo, uint64_t hi);
    double get_avg_coverage() const;

    //bool is_similar(double val1, double val2);


    // MEMBER VARIABLES
    // TODO ensure that there are never more than INT_MAX expressed regions per chromosome!
    std::vector<int> er_ids; //all expressed regions in a stitched_ER, er_id corresponds to index of averager.expressed_regions and is -1 for spliced regions
    // example: stitched_ER consists of er_ids 45, 46, 47, 49 == vector indices of expressed_regions
    double across_er_coverage; // avg (weighted) coverage of all exons that are part of the stitched ER so far
    std::vector<std::pair<unsigned int, double>> all_coverages; // stores a pair of er length (= weight) + normalized average coverage of the er
    // Genomic [start, end) of each entry in er_ids, aligned by index. For a
    // real ER it is the coverage-derived extent; for a spliced region (-1) it
    // is the splice junction's [donor, acceptor]. write_to_gtf snaps an exon
    // edge to the adjacent junction coordinate and falls back to the coverage
    // extent where an edge has no junction.
    std::vector<std::pair<uint64_t, uint64_t>> er_bounds;
    unsigned int total_length; // combined length of all ERs (excluding spliced regions)
    uint64_t start;
    uint64_t end;
    std::string chrom;
    // Strand of the splice junctions that produced this stitched ER. '.' for
    // unstranded inputs (the SJs in this chromosome were not partitioned by
    // strand) and '+' or '-' when the stitching pass restricted to one
    // strand. Matches the strand semantics on BedGraphRow and SJRow.
    char strand = '.';

    // overload output operator for SJRow
    friend std::ostream& operator<< (std::ostream& os, const StitchedER& stitched_er)
    {

        for (unsigned int i = 0; i < stitched_er.er_ids.size(); i++)
        {
            os << i << "\t" << "(" <<  stitched_er.all_coverages[i].first <<"," << stitched_er.all_coverages[i].second << ")" << std::endl;
        }
        return os << stitched_er.across_er_coverage << "\t" << stitched_er.start << "\t" << stitched_er.end << "\t" << stitched_er.total_length << std::endl;

    }
};

#endif //MLS_STITCHEDER_H
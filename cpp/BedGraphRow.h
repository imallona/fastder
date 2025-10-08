//
// Created by marti on 08/10/2025.
//

#ifndef FASTDER_BEDGRAPHROW_H
#define FASTDER_BEDGRAPHROW_H

#endif //FASTDER_BEDGRAPHROW_H

class BedGraphRow
{
    public:
        std::string chrom;
        int start;
        int end;
        double coverage;
        //double avg = 0;
        int total_reads = 0;
        int length = 0;
        // add optional values for average coverage, DER identifier

        //print BedGraphRow
        void print() const;

};
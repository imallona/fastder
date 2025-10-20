//
// Created by martinalavanya on 20.10.25.
//

#ifndef MLS_INTEGRATOR_H
#define MLS_INTEGRATOR_H
#include <unordered_map>

#endif //MLS_INTEGRATOR_H

#include <StitchedER.h>

class Integrator
{
    public:
    Integrator() = default;

    void stitch_up(std::unordered_map<unsigned int, unsigned int>& mm_sj_counts);


    // MEMBERS
    std::vector<StitchedER> stitched_ERs;

};
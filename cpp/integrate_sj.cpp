 // function that calculates relative match with a tolerance of +/- 5%
    bool is_similar(double val1, double val2){
      double tolerance_bottom = val1 * 0.95;
      double tolerance_top = val1 * 1.05;
      return val2 >= tolerance_bottom && val2 <= tolerance_top;
    }
    
    struct stitchedER
{
    std::vector<unsigned int> er_ids; //map to index of results
    double across_er_coverage; // avg (weighted) coverage of the stitched ER so far
    std::vector<pair<int, double>> all_coverages; // stores a pair of er length (= weight) + normalized average coverage of the er
    unsigned int total_reads = 0;
    unsigned int length = 0;
    // add optional values for average coverage, DER identifier

    // print BedGraphRow
    void print() const
    {
        std::cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << total_reads <<  "\t" << length << std::endl;
    }
    // later change to add weight / memory
    double get_avg_coverage(){
      double sum = 0;
      unsigned int total_length = 0;
      for (auto er : this->all_coverages){
        sum += er.first() * er.second(); // weighted sum of avg coverages across ERs
        total_length += er.second();
        
      }
      return sum / total_length;
    }

};
    
    
    
    
  for sample in samples:
    std::map<BedGraphRow> stitched_regions;
    std::vector<SJRow> sjs = mm_by_samples[sample.id] // all sjs that were discovered by STAR within one sample
    auto current_sj_id = sjs.begin(); //iterator over the vector of SJ ids
    std::vector<stitchedER> stitched_ers;
    StitchedER er1 = {0, {results[0].total_reads}, results[0].coverage}; //first StitchedER, currently consisting of 1 ER
    stitched_ers.push_back(er1); //CREATE STITCHED_ER FROM FIRST ER, vector has just one element so far
  
    // iterate over expressed regions, start at ER 2
    // check for each ER if it could belong to the most recent (stitched) ERs by checking if
    //
    for (unsigned int i = 1; i < results.size(); ++i){
      BedGraphRow& er = results[i];
      StitchedER& most_recent_er = stitched_ers.last(); //reference!
      // check if the current expressed region end coordinate somewhat matches a splice junction
      // check if current ER end coordinate is after the start of a SJ
      // do I use a tolerance region of +- 10 bp or something here?
      // check if end of sj and start of next er are in similar locations
       // check if coverage of current (stitched) er are similar -> likely a spliced region between two exons, stitch them together
        
      //TODO: make member function out of these three checks
      if (is_similar(most_recent_er.end, rr[current_sj_id].start) 
        && is_similar(result[i].start == rr[current_sj_id].end) 
        && is_similar(most_recent_er.across_er_coverage, result[i].coverage) {
        // check if next ER begins around end of the SJ
          most_recent_er.er_ids.push_back(i); // i is the index of the current element of results, i.e. the er_id
          most_recent_er.all_coverages.push_back(result[i].coverage);
          most_recent_er.across_er_coverage = most_recent_er.get_avg_coverage();

          // move to next sj id
          ++current_sj_id;
      }
      else if (!is_similar(most_recent_er.end, rr[current_sj_id].start) && most_recent_er.end > rr[current_sj_id].start){
        //TODO: checks if the next SJ might be a match with the current ER
        //TODO: should this be a while loops??
        append_er(++current_sj_id, er); // note that this increases current_sj_er as well!
      }
      
      // if neither match, there is no SJ right now, so we don't stitch
      // TODO do we instead create a new ER now?
      else {
        stitched_ers.push_back(er); //TODO modify this, maybe write constructor or member fct that does this cleanly
      }
      
      }
      // else: if really similar coverage, it could still be an exon! implement this for --mode adventure 
  }

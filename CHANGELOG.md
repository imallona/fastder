# Changelog

All notable changes to fastder are recorded in this file.

## [Unreleased]

### Fixed

- GTF output is now 1-based. fastder's internal coordinates are 0-based
  half-open (the BedGraph and BigWig convention), but the GTF writer emitted
  the start column without converting it, so every gene, transcript and exon
  start was one base too low. The start is now shifted by one on output. The
  end needs no shift, because a 0-based half-open end equals the 1-based
  inclusive end of the same interval.

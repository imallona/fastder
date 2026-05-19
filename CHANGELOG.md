# Changelog

All notable changes to fastder are recorded in this file.

## [0.1.0] - 2026-05-19

### Fixed

- GTF output is now 1-based. fastder's internal coordinates are 0-based
  half-open (the BedGraph and BigWig convention), but the GTF writer emitted
  the start column without converting it, so every gene, transcript and exon
  start was one base too low. The start is now shifted by one on output. The
  end needs no shift, because a 0-based half-open end equals the 1-based
  inclusive end of the same interval.
- RR splice junction coordinates are converted to 0-based on read. The RR
  file carries 1-based inclusive intron coordinates (the STAR SJ.out.tab
  convention), but fastder works in 0-based half-open coordinates like the
  coverage. The mismatch left junction-snapped exon ends one base too long.
  The RR parser now shifts the junction start by one on read, and rejects a
  junction start of 0, which is malformed under the 1-based convention,
  rather than letting the unsigned subtraction wrap.

### Changed

- Stitching is no longer gated on coverage similarity. The default
  `--coverage-tolerance` is raised to 1000, which makes the gate effectively
  off. A splice junction is the evidence that two regions belong to one
  transcript, and exon coverage varies widely within a transcript, so the
  previous default (0.8) wrongly rejected legitimate stitches. The flag is
  kept as a knob.
- Exon boundaries are snapped to splice junctions. When a stitched chain
  crosses a junction, the exon edges on either side are set to the junction's
  donor and acceptor coordinates instead of the coverage-derived expressed-
  region edges. Edges with no adjacent junction, namely the outer ends of a
  chain and both ends of a monoexonic region, keep the coverage extent.

## [Legacy]

Baseline forked from the upstream fastder repository, which carried no tagged
release. Versioned history begins at 0.1.0 above.

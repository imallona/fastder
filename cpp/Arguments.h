#ifndef FASTDER_ARGUMENTS_H
#define FASTDER_ARGUMENTS_H

#include <charconv>
#include <optional>
#include <string>
#include <string_view>

// Digits only, within the unsigned range.
inline std::optional<unsigned int> parse_count(std::string_view text)
{
    unsigned int value = 0;
    const char* end = text.data() + text.size();
    const auto parsed = std::from_chars(text.data(), end, value);
    if (text.empty() || parsed.ec != std::errc{} || parsed.ptr != end) return std::nullopt;
    return value;
}

// The last two options are named only when set, so earlier names stay valid.
inline std::string result_file_name(int position_tolerance, double min_coverage,
                                    double coverage_tolerance, int min_length,
                                    unsigned int min_junction_reads, bool no_stitch)
{
    std::string name = "FASTDER_RESULT_POS_TOL_" + std::to_string(position_tolerance)
        + "_MIN_COV_" + std::to_string(min_coverage)
        + "_COV_TOL_" + std::to_string(coverage_tolerance)
        + "_MIN_LENGTH_" + std::to_string(min_length);
    if (min_junction_reads > 0) name += "_MIN_JUNCTION_READS_" + std::to_string(min_junction_reads);
    if (no_stitch) name += "_NO_STITCH";
    return name + ".gtf";
}

#endif //FASTDER_ARGUMENTS_H

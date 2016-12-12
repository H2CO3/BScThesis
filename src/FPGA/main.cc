#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cstring>

#include <hls_stream.h>

#include "align.hh"


template<typename T>
void fill_stream(hls::stream<T> &stream, const T *begin, const T *end)
{
	while (begin < end) {
		stream.write(*begin++);
	}
}

std::vector<Dihedral> read_sequences(std::istream &instream)
{
    std::vector<Dihedral> vec;

    angle_type phi, psi;
    while (instream >> phi >> psi) {
        vec.push_back({ phi, psi });
    }

    return vec;
}

std::vector<index_type> read_lengths(std::istream &instream)
{
    std::vector<index_type> vec;

    std::string line;
    std::getline(instream, line);

    const char *ptr = line.data();
    const char *end = ptr + line.size();

    while (ptr < end) {
        char *next = nullptr;
        errno = 0;

        index_type length = std::strtol(ptr, &next, 10);

        if (length == 0 && next == ptr || errno != 0) {
            break;
        }

        vec.push_back(length);
        ptr = next;
    }

    return vec;
}



int main(int argc, char *argv[])
{
    // Parse arguments
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <scoring_offset> <gap_penalty>\n", argv[0]);
        return -1;
    }

    score_type scoring_offset = std::strtol(argv[1], nullptr, 10);
    score_type gap_penalty    = std::strtol(argv[2], nullptr, 10);

    // this is here so that the input stream can be changed easily later
    std::istream &instream = std::cin;

    // throw away explicit number of sequences
    {
        std::string line;
        std::getline(instream, line);
    }

    // Read lengths of sequences
    auto lengths = read_lengths(instream);

    // Read sequence data
    auto sequences = read_sequences(instream);

    // Perform alignment
    hls::stream<axi_out_score_type> out_scores;

    std::size_t ver_begin = 0;
    double elapsed_time = 0.0;

    using ull = unsigned long long;
    ull num_cells = 0;

    for (std::size_t i = 0; i < lengths.size() - 1; i++) {
        std::size_t ver_end = ver_begin + lengths[i];
        std::size_t hor_begin = ver_end;
        std::size_t hor_end = sequences.size();

        hls::stream<Dihedral>   stream_ver;
        hls::stream<Dihedral>   streams_hor;
        hls::stream<index_type> stream_sizes_hor;

        fill_stream(stream_ver,       &sequences[ver_begin], &sequences[ver_end]);
        fill_stream(streams_hor,      &sequences[hor_begin], &sequences[hor_end]);
        fill_stream(stream_sizes_hor, &lengths[i + 1],       &lengths[lengths.size()]);

        num_cells += (ull)(ver_end - ver_begin) * (ull) std::accumulate(&lengths[i + 1], &lengths[lengths.size()], 0ull);

        auto t_begin = std::chrono::steady_clock::now();

        align(
            stream_ver,
            stream_ver.size(),
            streams_hor,
            stream_sizes_hor,
            stream_sizes_hor.size(),
            scoring_offset,
            gap_penalty,
            out_scores
        );

        auto t_end = std::chrono::steady_clock::now();
        elapsed_time += std::chrono::duration<double>(t_end - t_begin).count();

        std::fprintf(stderr, "%lg... ", elapsed_time);
        std::fflush(stderr);

        ver_begin = ver_end;
    }

    // Dump performance counter to stderr
    std::fprintf(stderr, "\nElapsed time: %lg seconds\nNumber of cells: %llu\n", elapsed_time, num_cells);

    // Dump results
    std::size_t group_length = lengths.size() - 1;
    std::size_t group_index = 0;
    std::size_t score_index = 0;
    bool should_print_group_index = true;

	while (not out_scores.empty()) {
        if (should_print_group_index) {
            should_print_group_index = false;
            std::printf("#%zu.\t", group_index++);
        }

		std::printf(" %ld", static_cast<long>(out_scores.read().data));

        if (++score_index == group_length) {
            score_index = 0;
            --group_length;
            should_print_group_index = true;
            std::printf("\n");
        }
	}

	std::printf("\n");

    return 0;
}

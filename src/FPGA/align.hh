//
// align.hh
//
// Syntesizable (hopefully) part of the
// parallel Smith-Waterman algorithm
//
// Created on 06/03/2015
// by Arpad Goretity
//

#ifndef SWPARA_ALIGN_HH
#define SWPARA_ALIGN_HH

#include <cstdint>
#include <climits>
#include <hls_stream.h>
#include <ap_axi_sdata.h>

#include "util.hh"


#define MAX_SEQ_SIZE      512 // preferably a power of 2

// Horizontal and vertial size of sliding window
// TODO: update these (uncomment currently-commented values)
// once handling of shorter-than-full-length sequences is implemented.
#define WIN_COLS          16 // preferably a power of 2
#define WIN_ROWS          MAX_SEQ_SIZE // must be a power of 2
#define WIN_ROWS_MASK     (WIN_ROWS - 1)

static_assert(WIN_ROWS > 0 && (WIN_ROWS & (WIN_ROWS - 1)) == 0, "window height must be a power of two");

// Number of diagonal passes per sligind window
#define WIN_DIAGS         (WIN_ROWS + WIN_COLS - 1)

// Number of sliding windows horizontally and vertically
// TODO: update these (uncomment currently-commented values)
// once handling of shorter-than-full-length sequences is implemented.
#define WIN_COUNT_HOR     (MAX_SEQ_SIZE / WIN_COLS)
#define WIN_COUNT_VER     (MAX_SEQ_SIZE / WIN_ROWS)

// also try:
// * int17? (compiles only under synthesis!)
// * float
// ap_fixed<17, 1, 8, 8>?
// XXX: this MUST be a signed type.
typedef std::int16_t angle_type;

typedef std::uint16_t unsigned_angle_type;

static_assert(sizeof(angle_type) == sizeof(unsigned_angle_type), "signed and unsigned angle types must have the same size");

// TODO: what should this be? are 32 bits enough?
// XXX: this MUST be a signed type.
typedef std::int32_t score_type;

static_assert(sizeof(unsigned_angle_type) < sizeof(score_type), "signed score_type must be able to represent every value of unsigned_angle_type");

// This MUST be ap_axis<> if score_type is signed.
// This MUST be ap_axiu<> if score_type is unsigned.
typedef ap_axis<sizeof(score_type) * CHAR_BIT, 1, 1, 1> axi_out_score_type;

struct Dihedral {
	angle_type phi;
	angle_type psi;
};


extern "C" void align(
	hls::stream<Dihedral> &stream_ver,
	index_type stream_size_ver,
	hls::stream<Dihedral> &streams_hor,
	hls::stream<index_type> &stream_sizes_hor,
	seq_count_type num_streams_hor,
	score_type scoring_offset,
	score_type gap_penalty,
	hls::stream<axi_out_score_type> &out_scores
);

#endif // SWPARA_ALIGN_HH

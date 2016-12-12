//
// align.cc
//
// Synthesizable (hopefully) part of the
// parallel Smith-Waterman algorithm
//
// Created on 06/03/2015
// by Arpad Goretity
//

#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "util.hh"
#include "align.hh"


// The optimizer is not smart enough to realize that a loop has a constant
// iteration count when its begin and end iterators are known at compile-time.
// Consequently, using std::max(std::initializer_list) makes pipelining fail.
// Hence, I had to write this function in which the induction variable
// is compared to an integral constant, not an iterator...
template<typename T, std::size_t N>
static T array_max(std::array<T, N> arr)
{
#pragma HLS INLINE

	static_assert(N > 0, "cannot compute maximum of empty array");

	auto result = arr[0];

	for (std::size_t i = 1; i < arr.size(); i++) {
		if (arr[i] > result) {
			result = arr[i];
		}
	}

	return result;
}

static unsigned_angle_type angle_abs_diff(angle_type x, angle_type y)
{
#pragma HLS INLINE
	angle_type signed_diff = angle_type(unsigned_angle_type(x) - unsigned_angle_type(y));
	// An explicit implementation of absolute value is done here inline,
	// because std::abs<angle_type> compiles to crazy floating-point stuff.
	return signed_diff < 0 ? -signed_diff : +signed_diff;
}

static score_type dihedral_score(Dihedral angle1, Dihedral angle2, score_type offset)
{
#pragma HLS INLINE
	score_type dphi = angle_abs_diff(angle1.phi, angle2.phi);
	score_type dpsi = angle_abs_diff(angle1.psi, angle2.psi);
	return offset - (dphi * dphi + dpsi * dpsi);
}

score_type align_one(
	hls::stream<Dihedral> &stream_ver,
	index_type stream_size_ver,
	hls::stream<Dihedral> &stream_hor,
	index_type stream_size_hor,
	score_type scoring_offset,
	score_type gap_penalty,
	bool should_read_ver_stream
)
{
#pragma HLS INLINE

	assert(WIN_COLS >= 2 && "window columns >= 2 required for horizontal propagation to work correctly "
				"(if you need less than 2 columns, add a temporary horizontal propagation buffer)");


	// in the following declarations, -1 means 'Uninitialized'
#ifdef __SYNTHESIS__
	static Dihedral seq_ver[WIN_ROWS];
	static Dihedral seq_hor[WIN_COLS];

#pragma HLS reset variable=seq_ver off
#pragma HLS reset variable=seq_hor off

	score_type hor_prop_buf[WIN_ROWS]; // = { 0 };

	score_type diag_buf_old[WIN_COLS];
	score_type diag_buf_new[WIN_COLS];
#else
	static std::vector<Dihedral> seq_ver(WIN_ROWS, { -1, -1 });
	static std::vector<Dihedral> seq_hor(WIN_COLS, { -1, -1 });

	std::vector<score_type> hor_prop_buf(WIN_ROWS, -1);

	// 1000 is an arbitrarily big positive pseudo-"garbage" value that is
	// used for checking whether boundary conditions are implemented correctly
	// so that huge leftover values don't mess up computation of the maximum.
	std::vector<score_type> diag_buf_old(WIN_COLS, 1000);
	std::vector<score_type> diag_buf_new(WIN_COLS, 1000);
#endif

	assert(std::end(seq_ver) - std::begin(seq_ver) == WIN_ROWS);
	assert(std::end(seq_hor) - std::begin(seq_hor) == WIN_COLS);
	assert(WIN_ROWS > WIN_COLS && "window must be higher than wide");

	hls_debug("stream_ver size: %zu\n", stream_ver.size());
	hls_debug("stream_hor size: %zu\n", stream_hor.size());

	// h: index of horizontal (rightward) sliding window
	// i: index of the diagonal within the current window (transformed index)
	// j: The index of the current cell within the diagonal it lies in (transformed index)
	// r: index of the row of the current cell within the current window (non-transformed index)
	// c: index of the column of the current cell within the current window (non-transformed index)
	index_type h, i, j, r, c;

	// These are set to the appropriate fraction of the window width and height,
	// respectively, when the sequence size is not an integer multiple thereof.
	index_type max_valid_row;
	index_type max_valid_col;

	// +------+------+
	// |(0, 0)|(0, 1)|
	// +------+------+
	// |(1, 0)|(1, 1)|
	// +------+------+
	//
	// This 2 x 2 lookahead window is read and written all the time,
	// so we want it to be synthesized as simple registers, not as block BRAM.
	//
	// The negative powers of two serve merely as distinct-from-negative-one
	// indicators of "uninitialized value". This window is updated in a manner
	// so that it doesn't need to be initialized explicitly.
	score_type lah_buf[2][2]; // = { { -128, -256 }, { -512, -1024 } };
#pragma HLS ARRAY_PARTITION variable=lah_buf complete dim=0

	// we always read the next element of the diagonal buffers into
	// registers so that they can be used as many times as necessary.
	score_type diag_buf_old_next_cell;
	score_type diag_buf_new_next_cell;

	// We should not read the horizontal propagation buffer twice in an iteration,
	// hence we created this small shift register.
	//
	// +---+
	// | 0 | i - 2   |
	// +---+         |  sliding downward
	// | 1 | i - 1   v
	// +---+
	score_type hor_prop_buf_next_cells[2]; // = { -2048, -4096 };
#pragma HLS ARRAY_PARTITION variable=hor_prop_buf_next_cells complete dim=0

	// score is always non-negative -> this is OK
	score_type max_score = 0;
	score_type cur_score;

	// These registers are used for reading into seq_hor and ver_hor RAMs
	Dihedral seq_hor_read_reg;
	Dihedral seq_ver_read_reg;

	// ...and these two complement the ones above: they serve as temporaries
	// for the actual score computation.
	Dihedral seq_hor_comp_reg;
	Dihedral seq_ver_comp_reg;

	// Holds the original size of the vertical input stream.
	const index_type stream_size_ver_orig = stream_size_ver;

	// this indicates whether the row index is within bounds (0 <= r < WIN_ROWS),
	// so it effectively serves as a 'Write Enable' signal for hor_prop_buf and seq_ver.
	bool in_bounds;

hor_window_loop:
	for (h = 0; h < WIN_COUNT_HOR; h++) {

	diag_loop:
		for (i = 0; i < WIN_DIAGS; i++) {

		col_loop:
			for (j = 0; j < WIN_COLS; j++) {
#pragma HLS DEPENDENCE variable=seq_hor false
#pragma HLS DEPENDENCE variable=seq_ver false
#pragma HLS DEPENDENCE variable=hor_prop_buf false
#pragma HLS PIPELINE II=1 rewind

				// compute non-transformed indices from transformed ones
				r = i - j;
				c = j;

				// read horizontal propagation buffer at the beginning of each diagonal,
				// and update its temporary shift register
				if (j == 0 && i < WIN_ROWS) {
					hor_prop_buf_next_cells[0] = 0 < i ? hor_prop_buf_next_cells[1] : 0;
					hor_prop_buf_next_cells[1] = 0 < h ? hor_prop_buf[i] : 0;
					hls_debug("hor_prop_buf_next_cells = [%d, %d]\n", hor_prop_buf_next_cells[0], hor_prop_buf_next_cells[1]);
				}

				// Update temporary registers
				if (j == 0) {
					lah_buf[0][0] = i < 1 ? 0 : hor_prop_buf_next_cells[0];
					lah_buf[1][0] = i < 0 ? 0 : hor_prop_buf_next_cells[1];
				} else {
					lah_buf[0][0] = lah_buf[0][1];
					lah_buf[1][0] = lah_buf[1][1];
				}

				// Read ahead, respecting boundary conditions.
				// Out-of-bounds values should be 0, so that 'partial'
				// windows (potentially occurring at the end of sequences
				// when the sequence length is not an integer multiple of
				// the window size) will have correct values even on the
				// boundaries (since not all invalid elements of a partial
				// window are out of bounds!)
				diag_buf_old_next_cell = j <= i - 2 ? diag_buf_old[j] : 0;
				diag_buf_new_next_cell = j <= i - 1 ? diag_buf_new[j] : 0;

				lah_buf[0][1] = i < 2 ? 0 : diag_buf_old_next_cell;
				lah_buf[1][1] = i < 1 ? 0 : diag_buf_new_next_cell;

				// if the cell coordinates are OOB, don't try to compute them.
				// 'c' should always be within bounds, because it's equal to j.
				in_bounds = 0 <= r && r < WIN_ROWS /* && 0 <= c && c < WIN_COLS */;
				assert(0 <= c && c < WIN_COLS);

				// Set row and column boundaries in order to handle partial windows
				if (i == 0 && j == 0) {
					if (h == 0) {
						max_valid_row = WIN_ROWS;
					}

					max_valid_col = WIN_COLS;
				}

				// Read in vertical sequence buffer if necessary,
				// i.e. at the beginning of every horizontal row of windows,
				// indicated by the horizontal window index being reset to 0.
				if (h == 0 && j == 0) {
					// Only read if there are any elements left in the stream.
					// Otherwise, register the end-of-stream state.
					// hls::stream::empty() doesn't do in hardware what one might think it should do...
					// Hence, an explicit check-size-and-decrement is performed.
					// if (stream_ver.empty()) {
					if (stream_size_ver == 0) {
						if (r < max_valid_row) {
							max_valid_row = r;
							hls_debug("setting max_valid_row = %td\n", std::ptrdiff_t(max_valid_row));
						}
					} else {
						if (in_bounds) {
							if (should_read_ver_stream) {
								seq_ver_read_reg = stream_ver.read();
								seq_ver[r] = seq_ver_read_reg;
							}

							stream_size_ver--;
						}

						hls_debug(
							"reading seq_ver[%td] = (%d, %d)\n",
							std::ptrdiff_t(r),
							seq_ver[size_type(r) & WIN_ROWS_MASK].phi,
							seq_ver[size_type(r) & WIN_ROWS_MASK].psi
						);
					}
				}

				// Read in horizontal sequence buffer if necessary,
				// i.e. at the beginning of each window
				if (r == 0) {
					// if (stream_hor.empty()) {
					if (stream_size_hor == 0) {
						if (c < max_valid_col) {
							max_valid_col = c;
							hls_debug("setting max_valid_col = %td\n", std::ptrdiff_t(max_valid_col));
						}
					} else {
						seq_hor_read_reg = stream_hor.read();
						seq_hor[c] = seq_hor_read_reg;
						stream_size_hor--;

						hls_debug("reading seq_hor[%td] = (%d, %d)\n", std::ptrdiff_t(c), seq_hor[c].phi, seq_hor[c].psi);
					}
				}

				// Compute score
				// (as long as we are within bounds of valid sequence data)
				// We can just play the "compute invalid but ignore" game,
				// since data dependencies are monotonic: each cell depends
				// on previous cells (those with smaller indices) only, so
				// calculating garbage values doesn't affect the correctness
				// of valid, within-bounds cells.
				if (r == 0) {
					seq_hor_comp_reg = seq_hor_read_reg;
				} else {
					seq_hor_comp_reg = seq_hor[c];
				}

				if (h == 0 && j == 0 && should_read_ver_stream) {
					seq_ver_comp_reg = seq_ver_read_reg;
				} else {
					seq_ver_comp_reg = seq_ver[size_type(r) & WIN_ROWS_MASK];
				}

				cur_score = array_max<score_type, 4>({
					lah_buf[0][0] /* diag neighbor */ + dihedral_score(seq_ver_comp_reg, seq_hor_comp_reg, scoring_offset),
					lah_buf[1][0] /* left neighbor */ + gap_penalty,
					lah_buf[1][1] /* top  neighbor */ + gap_penalty,
					score_type(0) /* align locally */
				});

				// accumulate maximum if it's valid
				if (in_bounds && r < max_valid_row && c < max_valid_col) {
					if (cur_score > max_score) {
						max_score = cur_score;
					}

					hls_debug("r=%td  c=%td\n", std::ptrdiff_t(r), std::ptrdiff_t(c));
					hls_debug("%3d %3d\n%3d %3d\n", lah_buf[0][0], lah_buf[0][1], lah_buf[1][0], lah_buf[1][1]);
					hls_debug("    score = %d\n", cur_score);
				}

				// Shift values in diagonal buffers
				diag_buf_old[j] = diag_buf_new_next_cell;
				diag_buf_new[j] = cur_score;

				// propagate cells in the last column of each row rightwards
				if (c == WIN_COLS - 1) {
					hls_debug("    rightward-propagating end of row[%td] = %d\n", std::ptrdiff_t(r), cur_score);

					if (in_bounds) {
						hor_prop_buf[r] = cur_score;
					}
				}

				hls_debug("\n");
			}

			if (i == stream_size_ver_orig + WIN_COLS - 1) {
				break;
			}
		}

		if (stream_size_hor == 0) {
			break;
		}
	}

	return max_score;
}

void align(
	hls::stream<Dihedral> &stream_ver,
	index_type stream_size_ver,
	hls::stream<Dihedral> &streams_hor,
	hls::stream<index_type> &stream_sizes_hor,
	seq_count_type num_streams_hor,
	score_type scoring_offset,
	score_type gap_penalty,
	hls::stream<axi_out_score_type> &out_scores
)
{
#pragma HLS INTERFACE s_axilite port=stream_size_ver
#pragma HLS INTERFACE s_axilite port=num_streams_hor
#pragma HLS INTERFACE s_axilite port=scoring_offset
#pragma HLS INTERFACE s_axilite port=gap_penalty

#pragma HLS INTERFACE axis port=stream_ver
#pragma HLS DATA_PACK variable=stream_ver

#pragma HLS INTERFACE axis port=streams_hor
#pragma HLS DATA_PACK variable=streams_hor

#pragma HLS INTERFACE axis port=stream_sizes_hor

#pragma HLS INTERFACE axis port=out_scores

#pragma HLS INTERFACE s_axilite port=return

	bool should_read_ver_stream = true;

align_loop:
	for (seq_count_type i = 0; i < num_streams_hor; i++) {
		index_type stream_size_hor = stream_sizes_hor.read();

		// Compute score
		score_type score = align_one(
			stream_ver,
			stream_size_ver,
			streams_hor,
			stream_size_hor,
			scoring_offset,
			gap_penalty,
			should_read_ver_stream
		);

		// Convert raw numeric value to AXI streamable type with side-band signals.
		// Set the TLAST bit so that the DMA knows when to flush a potential partial burst.
		// Set keep and strobe signals to all 1's (-1 in 2's complement)
		auto axi_score = axi_out_score_type{};

		axi_score.data = score;
		axi_score.keep = -1;
		axi_score.strb = -1;
		axi_score.last = i == num_streams_hor - 1;

		out_scores.write(axi_score);

		// only read vertical stream upon alignment of first horizontal sequence
		should_read_ver_stream = false;
	}
}

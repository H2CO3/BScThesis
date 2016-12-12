/*
 * align_fpga.h
 *
 *  Created on: Oct 29, 2016
 *      Author: h2co3
 *
 * Driver for alignment hardware implemented on FPGA
 */

#ifndef ALIGN_FPGA_H_
#define ALIGN_FPGA_H_

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

#include "platform.h"
#include "xalign.h"
#include "xaxidma.h"
#include "ff.h"


// Log scores only for debugging purposes and small databases.
// If the system is used with a large database, where the number of sequences
// is in the ten thousands or hundred thousands, then outputting a quadratic
// amount of scores is not feasible through the slow USART port.
#define USART_LOG_SCORES 1


typedef  int16_t index_type;
typedef uint16_t size_type;

typedef  int16_t angle_type;
typedef  int32_t score_type;

typedef uint32_t seq_count_type;

typedef struct Dihedral {
	angle_type phi;
	angle_type psi;
} Dihedral;

typedef struct Sequences {
	Dihedral *buffer;             // Raw sequence data, contiguously (owning pointer)
	index_type *sequence_lengths; // Number of dihedrals in each sequence (owning pointer)
	seq_count_type num_sequences; // Number of sequences
} Sequences;

typedef struct AlignSystem {
	XAxiDma ver_axidma;
	XAxiDma hor_axidma;
	XAxiDma hor_sizes_axidma;
	XAxiDma out_scores_axidma;
	XAlign align;
} AlignSystem;


void init_align_system(AlignSystem *align_sys, score_type scoring_offset, score_type gap_penalty);

size_t total_seq_len(const index_type *seq_lens, seq_count_type num_seqs);

bool run_align(
	AlignSystem *align_sys,
	Sequences *seqs,
	FIL *out_file,
	double *elapsed_time
);

#endif /* ALIGN_FPGA_H_ */

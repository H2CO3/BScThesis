/*
 * align_fpga.c
 *
 *  Created on: Oct 29, 2016
 *      Author: h2co3
 *
 * Driver for alignment hardware implemented on FPGA
 */

#include <inttypes.h>

#include "xtime_l.h"
#include "align_fpga.h"
#include "seq_file.h"


#define VER_AXIDMA_ID        XPAR_AXI_DMA_1_DEVICE_ID
#define HOR_AXIDMA_ID        XPAR_AXI_DMA_0_DEVICE_ID
#define HOR_SIZES_AXIDMA_ID  XPAR_AXI_DMA_2_DEVICE_ID
#define OUT_SCORES_AXIDMA_ID XPAR_AXI_DMA_3_DEVICE_ID
#define ALIGN_ID             XPAR_ALIGN_0_DEVICE_ID


static void init_align(XAlign *instance, u32 device_id)
{
	XAlign_Config *config = XAlign_LookupConfig(device_id);
	XAlign_CfgInitialize(instance, config);
}

static void init_axidma(XAxiDma *instance, u32 device_id)
{
	XAxiDma_Config *config = XAxiDma_LookupConfig(device_id);
	XAxiDma_CfgInitialize(instance, config);
	XAxiDma_IntrDisable(instance, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DEVICE_TO_DMA);
	XAxiDma_IntrDisable(instance, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DMA_TO_DEVICE);
}

void init_align_system(AlignSystem *align_sys, score_type scoring_offset, score_type gap_penalty)
{
	init_align (&align_sys->align,             ALIGN_ID);
	init_axidma(&align_sys->ver_axidma,        VER_AXIDMA_ID);
	init_axidma(&align_sys->hor_axidma,        HOR_AXIDMA_ID);
	init_axidma(&align_sys->hor_sizes_axidma,  HOR_SIZES_AXIDMA_ID);
	init_axidma(&align_sys->out_scores_axidma, OUT_SCORES_AXIDMA_ID);

	XAlign_Set_scoring_offset(&align_sys->align, scoring_offset);
	XAlign_Set_gap_penalty   (&align_sys->align, gap_penalty);
}

static void flush_cache(const void *addr, size_t elem_size, size_t num_elems)
{
	Xil_DCacheFlushRange((uintptr_t)addr, elem_size * num_elems);
}

static void invalidate_cache(void *addr, size_t elem_size, size_t num_elems)
{
	Xil_DCacheInvalidateRange((uintptr_t)addr, elem_size * num_elems);
}

static u32 axidma_write(XAxiDma *axidma, const void *buf, size_t elem_size, size_t num_elems)
{
	return XAxiDma_SimpleTransfer(
		axidma,
		(uintptr_t)buf,
		elem_size * num_elems,
		XAXIDMA_DMA_TO_DEVICE
	);
}

static u32 axidma_read(XAxiDma *axidma, void *buf, size_t elem_size, size_t num_elems)
{
	return XAxiDma_SimpleTransfer(
		axidma,
		(uintptr_t)buf,
		elem_size * num_elems,
		XAXIDMA_DEVICE_TO_DMA
	);
}

static bool axidma_busy_writing(XAxiDma *axidma)
{
	return XAxiDma_Busy(axidma, XAXIDMA_DMA_TO_DEVICE);
}

static bool axidma_busy_reading(XAxiDma *axidma)
{
	return XAxiDma_Busy(axidma, XAXIDMA_DEVICE_TO_DMA);
}

static bool align_busy(XAlign *align)
{
	return XAlign_IsDone(align) == false;
}

static XTime get_time(void)
{
	XTime t = 0;
	XTime_GetTime(&t);
	return t;
}

size_t total_seq_len(const index_type *seq_lens, seq_count_type num_seqs)
{
	size_t len = 0;

	for (seq_count_type i = 0; i < num_seqs; i++) {
		len += seq_lens[i];
	}

	return len;
}

bool run_align(
	AlignSystem *align_sys,
	Sequences *seqs,
	FIL *out_file,
	double *elapsed_time
)
{
	// Initialize timing info needed for benchmarking
	XTime total_delta_t = 0;
	*elapsed_time = 0.0; // just in case there's an error or some other early return

	// Compute total length of sequences
	size_t seq_len = total_seq_len(seqs->sequence_lengths, seqs->num_sequences);

	// Needed for computing the padding at the end of the output file
	size_t total_bytes_written = 0;

	// Write number of sequences to output file
	CHK_FOP(f_write_chk(out_file, &seqs->num_sequences, sizeof seqs->num_sequences));
	total_bytes_written += sizeof seqs->num_sequences;

	// Flush/invalidate caches under buffers explicitly for coherence
	flush_cache(seqs->buffer,           sizeof seqs->buffer[0],           seq_len);
	flush_cache(seqs->sequence_lengths, sizeof seqs->sequence_lengths[0], seqs->num_sequences);

	// Allocate memory for results
	score_type *out_scores = malloc(seqs->num_sequences * sizeof out_scores[0]);

	// Initialize indices, pointers and lengths that determine the range
	// of the vertical (resident) sequence and that of the horizontal (volatile) sequences
	size_t num_seqs_hor = seqs->num_sequences;
	size_t len_hor = seq_len;
	Dihedral *seq_ver = seqs->buffer;

	// Compare every sequence to every other one.
	// First, the vertical sequence is taken to be sequence #0, and it's compared
	// to sequences #1...n-1, which act as horizontal sequences.
	// Then, pointers, indices and lengths are updated so that sequence #1 becomes
	// the vertical one and it's being compared to seqs #2...n-1.
	// This repeats until the vertical sequence would be seq. #n-1, which is when the computation ends.
	for (seq_count_type i = 0; i < seqs->num_sequences - 1; i++) {
		index_type len_ver = seqs->sequence_lengths[i];
		Dihedral *seqs_hor = seq_ver + len_ver;

		len_hor -= len_ver;
		num_seqs_hor--;

		// invalidate part of cache where scores will be written
		invalidate_cache(out_scores, sizeof out_scores[0], num_seqs_hor);

		// Set stream lengths
		XAlign_Set_stream_size_ver(&align_sys->align, len_ver);
		XAlign_Set_num_streams_hor(&align_sys->align, num_seqs_hor);

		// Actually send the data
		u32 status = XST_SUCCESS;

#define CHECK(str) do { if (status != XST_SUCCESS) { printf("%s: status = %lu\r\n", str, status); free(out_scores); return false; } } while (0)

		status = axidma_write(
			&align_sys->ver_axidma,
			seq_ver,
			sizeof seq_ver[0],
			len_ver
		);
		CHECK("vertical sequence data");

		status = axidma_write(
			&align_sys->hor_axidma,
			seqs_hor,
			sizeof seqs_hor[0],
			len_hor
		);
		CHECK("horizontal sequence data");

		status = axidma_write(
			&align_sys->hor_sizes_axidma,
			&seqs->sequence_lengths[i + 1],
			sizeof seqs->sequence_lengths[i + 1],
			num_seqs_hor
		);
		CHECK("horizontal sequence lengths");

		status = axidma_read(
			&align_sys->out_scores_axidma,
			out_scores,
			sizeof out_scores[0],
			num_seqs_hor
		);
		CHECK("out scores");

#undef CHECK

		// Start alignment block
		XAlign_Start(&align_sys->align);

		// Wait for them to finish using polling.
		// Measure the elapsed time.
		XTime t_begin = get_time();

		while (
		     axidma_busy_writing(&align_sys->ver_axidma)
		  || axidma_busy_writing(&align_sys->hor_axidma)
		  || axidma_busy_writing(&align_sys->hor_sizes_axidma)
		  || axidma_busy_reading(&align_sys->out_scores_axidma)
		  || align_busy(&align_sys->align)
		) {
			// NOP
		}

		XTime t_end = get_time();
		total_delta_t += t_end - t_begin;

		// Dump scores via USART if necessary
#if USART_LOG_SCORES
		printf("#%" PRIu32 ".\t", i);

		for (size_t j = 0; j < num_seqs_hor; j++) {
			printf(" %" PRIi32, out_scores[j]);
		}

		printf("\r\n");
#endif

		// Write scores to file
		size_t score_bufsize = num_seqs_hor * sizeof out_scores[0];
		CHK_FOP(f_write_chk(out_file, out_scores, score_bufsize));
		total_bytes_written += score_bufsize;

		seq_ver += len_ver;
	}

	// Must pad file by rounding up to the maximal sector size,
	// otherwise nothing is written to the file whatsoever
	CHK_FOP(pad_file(out_file, total_bytes_written));
	CHK_FOP(f_sync(out_file));

	// Write performance info to out parameter
	*elapsed_time = total_delta_t * 1.0 / COUNTS_PER_SECOND;

	free(out_scores);
	return true;
}

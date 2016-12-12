/*
 * seq_file.c
 *
 *  Created on: Oct 29, 2016
 *      Author: h2co3
 *
 * Driver for reading sequences from / writing results to the SD Card
 */

#include "seq_file.h"


const char *const fresult_strings[] = {
	"Succeeded",
	"A hard error occurred in the low-level disk I/O layer",
	"Internal error / Assertion failed",
	"The physical drive does not work",
	"Could not find the file",
	"Could not find the path",
	"The path name format is invalid",
	"Access denied due to prohibited access or directory full",
	"Access denied due to prohibited access",
	"The file/directory object is invalid",
	"The physical drive is write-protected",
	"The logical drive number is invalid",
	"The volume has no work area",
	"There is no valid FAT volume",
	"The function f_mkfs() aborted due to a parameter error",
	"Could not get a grant to access the volume within the defined period",
	"The operation is rejected according to the file sharing policy",
	"LFN working buffer could not be allocated",
	"Number of open files greater than _FS_SHARE",
	"Invalid parameter",
};


FRESULT mount_fat_fs(FATFS *outfs)
{
	return f_mount(outfs, "0:/", 1);
}

FRESULT unmount_fat_fs(void)
{
	return f_mount(NULL, "0:", 0);
}

FRESULT f_read_chk(FIL *file, void *buf, size_t size)
{
	UINT bytes_read = 0;
	FRESULT result = f_read(file, buf, size, &bytes_read);

	if (result != FR_OK) {
		return result;
	}

	return bytes_read == size ? FR_OK : FR_INT_ERR;
}

FRESULT f_write_chk(FIL *file, const void *buf, size_t size)
{
	UINT bytes_written = 0;
	FRESULT result = f_write(file, buf, size, &bytes_written);

	if (result != FR_OK) {
		return result;
	}

	return bytes_written == size ? FR_OK : FR_INT_ERR;
}

FRESULT pad_file(FIL *file, size_t total_bytes_written)
{
	size_t sector_size = sizeof file->fs->win;
	size_t last_partial_chunk_size = total_bytes_written % sector_size;

	if (last_partial_chunk_size) {
		size_t trailing_size = sector_size - last_partial_chunk_size;
		unsigned char zeros[trailing_size];
		memset(zeros, 0, sizeof zeros);
		return f_write_chk(file, zeros, sizeof zeros);
	}

	return FR_OK;
}

FRESULT read_sequences_from_file(FIL *file, Sequences *seqs)
{
	FRESULT fresult = FR_OK;

	// First, read the number of sequences in the file
	seq_count_type num_seqs = 0;

	if ((fresult = f_read_chk(file, &num_seqs, sizeof num_seqs)) != FR_OK) {
		return fresult;
	}

	// Then read the length of each sequence
	index_type *seq_lens = NULL;
	size_t seq_lens_bufsize = num_seqs * sizeof seq_lens[0];
	seq_lens = malloc(seq_lens_bufsize);

	if ((fresult = f_read_chk(file, seq_lens, seq_lens_bufsize)) != FR_OK) {
		free(seq_lens);
		return fresult;
	}

	// Finally, read actual sequence data
	Dihedral *buf = NULL;
	size_t seq_bufsize = total_seq_len(seq_lens, num_seqs) * sizeof buf[0];
	buf = malloc(seq_bufsize);

	if ((fresult = f_read_chk(file, buf, seq_bufsize)) != FR_OK) {
		free(seq_lens);
		free(buf);
		return fresult;
	}

	// Populate out parameter
	seqs->buffer = buf;
	seqs->sequence_lengths = seq_lens;
	seqs->num_sequences = num_seqs;

	return FR_OK;
}

void free_sequences(Sequences *seqs)
{
	free(seqs->buffer);
	free(seqs->sequence_lengths);
}

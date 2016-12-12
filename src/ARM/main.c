/*
 * main.c: driver for the parallel Smith-Waterman alignment block
 *
 * This application configures UART 16550 to baud rate 9600.
 * PS7 UART (Zynq) is not initialized by this application, since
 * bootrom/bsp configures it to baud rate 115200
 *
 * ------------------------------------------------
 * | UART TYPE   BAUD RATE                        |
 * ------------------------------------------------
 *   uartns550   9600
 *   uartlite    Configurable only in HW design
 *   ps7_uart    115200 (configured by bootrom/bsp)
 */

#include "align_fpga.h"
#include "seq_file.h"


// Files to read sequences from and write results to
#define INPUT_FILENAME  "INPUT.BIN"
#define OUTPUT_FILENAME "OUTPUT.BIN"

// Parameters for the algorithm
#define SCORING_OFFSET  65536
#define GAP_PENALTY     (-4000)


int main()
{
	// Initialize black box stuff
	init_platform();
	printf("\r\n\r\n*** Initialized system!\r\n");

	// Mount FAT filesystem on SD Card
	FATFS fs;
	CHK_FOP(mount_fat_fs(&fs));
	printf("*** Mounted SD Card\r\n");

	// Open input file and read its contents
	FIL infile;
	CHK_FOP(f_open(&infile, INPUT_FILENAME, FA_READ));
	printf("*** Opened input file '%s'\r\n", INPUT_FILENAME);

	// Read sequence data into memory
	Sequences seqs;
	CHK_FOP(read_sequences_from_file(&infile, &seqs));
	printf("*** Read sequence data from file\r\n");

	// Close input file
	CHK_FOP(f_close(&infile));

	// Open output file
	FIL outfile;
	CHK_FOP(f_open(&outfile, OUTPUT_FILENAME, FA_WRITE | FA_OPEN_ALWAYS));

	// Initialize FPGA hardware
	AlignSystem align_sys;
	init_align_system(&align_sys, SCORING_OFFSET, GAP_PENALTY);
	printf("*** Initialized alignment hardware\r\n");

	// Compute results
	printf("*** Performing computations\r\n");

	double dt = 0.0;
	bool success = run_align(
		&align_sys,
		&seqs,
		&outfile,
		&dt
	);

	printf("*** %s! Elapsed Time: %lg seconds\r\n", success ? "Success" : "Failure", dt);

	// Finish writing results and close output file
	CHK_FOP(f_close(&outfile));

	printf("*** Results written to file '%s'\r\n", OUTPUT_FILENAME);

	// Clean up and exit
	free_sequences(&seqs);
	CHK_FOP(unmount_fat_fs());
	cleanup_platform();

	printf("\r\n");

	return 0;
}

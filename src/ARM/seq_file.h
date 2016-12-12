/*
 * seq_file.h
 *
 *  Created on: Oct 29, 2016
 *      Author: H2CO3
 *
 * Driver for reading sequences from / writing results to the SD Card
 */

#ifndef SEQ_FILE_H_
#define SEQ_FILE_H_

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

#include "align_fpga.h"
#include "ff.h"


#define CHK_FOP(fop)                                                        \
  do {                                                                      \
    FRESULT _fop_result_ = (fop);                                           \
    if (_fop_result_ != FR_OK) {                                            \
      printf("*** error: %s: %s\r\n", #fop, fresult_strings[_fop_result_]); \
      abort();                                                              \
    }                                                                       \
  } while (0)


extern const char *const fresult_strings[];


FRESULT mount_fat_fs(FATFS *outfs);
FRESULT unmount_fat_fs(void);

FRESULT f_read_chk(FIL *file, void *buf, size_t size);
FRESULT f_write_chk(FIL *file, const void *buf, size_t size);

FRESULT pad_file(FIL *file, size_t total_bytes_written);

FRESULT read_sequences_from_file(FIL *file, Sequences *seqs);

void free_sequences(Sequences *seqs);

#endif /* SEQ_FILE_H_ */

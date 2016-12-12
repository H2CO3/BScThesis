#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


typedef  int16_t angle_type;
typedef  int16_t index_type;
typedef  int32_t score_type;
typedef uint32_t seq_count_type;

typedef struct Dihedral {
    angle_type phi;
    angle_type psi;
} Dihedral;


static void generate(FILE *file);
static void dump_seq(FILE *file);
static void dump_score(FILE *file);


int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s [genseq | dumpseq | dumpscore] <file>\n", argv[0]);
        return -1;
    }

    if (strcmp(argv[1], "genseq") == 0) {
        FILE *file = fopen(argv[2], "wb");
        generate(file);
        fclose(file);
    } else if (strcmp(argv[1], "dumpseq") == 0) {
        FILE *file = fopen(argv[2], "rb");
        dump_seq(file);
        fclose(file);
    } else if (strcmp(argv[1], "dumpscore") == 0) {
        FILE *file = fopen(argv[2], "rb");
        dump_score(file);
        fclose(file);
    } else {
        fprintf(stderr, "Unknown command: %s\n", argv[1]);
        return -1;
    }

    return 0;
}

static void generate(FILE *file)
{
    seq_count_type num_seqs = arc4random_uniform(999) + 2; // 2...1000 sequences
    index_type seq_lens[num_seqs];

    for (seq_count_type i = 0; i < num_seqs; i++) {
        seq_lens[i] = arc4random_uniform(512); // sequences of maximal length 511
    }

    fwrite(&num_seqs, sizeof num_seqs, 1, file);
    fwrite(seq_lens, sizeof seq_lens[0], num_seqs, file);

    for (seq_count_type i = 0; i < num_seqs; i++) {
        for (index_type j = 0; j < seq_lens[i]; j++) {
            uint32_t tmp = arc4random();
            angle_type low  = (tmp >>  0) & 0xffff; // -32768...32767
            angle_type high = (tmp >> 16) & 0xffff; // -32768...32767
            Dihedral d = { low, high };
            fwrite(&d, sizeof d, 1, file);
        }
    }
}

static void dump_seq(FILE *file)
{
    seq_count_type num_seqs = 0;
    fread(&num_seqs, sizeof num_seqs, 1, file);

    index_type seq_lens[num_seqs];
    fread(seq_lens, sizeof seq_lens[0], num_seqs, file);

    printf("%lu\n", (unsigned long) num_seqs);

    for (seq_count_type i = 0; i < num_seqs; i++) {
        printf("%ld  ", (long) seq_lens[i]);
    }

    printf("\n");

    for (seq_count_type i = 0; i < num_seqs; i++) {
        for (index_type j = 0; j < seq_lens[i]; j++) {
            Dihedral d;
            fread(&d, sizeof d, 1, file);
            printf("%d %d      ", (int) d.phi, (int) d.psi);
        }

        printf("\n");
    }
}

static void dump_score(FILE *file)
{
    seq_count_type num_seqs = 0;
    fread(&num_seqs, sizeof num_seqs, 1, file);

    for (seq_count_type i = 0; i < num_seqs - 1; i++) {
        printf("#%lu.\t", (unsigned long) i);

        for (seq_count_type j = i + 1; j < num_seqs; j++) {
            score_type next_score;
            fread(&next_score, sizeof next_score, 1, file);
            printf(" %ld", (long) next_score);
        }

        printf("\n");
    }

    printf("\n");
}

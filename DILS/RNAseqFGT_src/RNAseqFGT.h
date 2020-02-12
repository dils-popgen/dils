#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE 100000


void get_IDs(char *string, char *ID_locus, char *ID_species, char *ID_individual, char *ID_allele);
void run_analysis(FILE *fasta_file, off_t pos, FILE * out_file, char *ID_locus, char *ID_species, int verbose);
int check_contig_sequences(FILE *fasta_file, off_t pos, off_t *last_pos, char *ID_locus, int *Lseq);
void check_sequence(char *string, char *ID_locus, char *ID_species, char *ID_individual, char *ID_allele);
char **load_sequence(FILE *fasta_file, off_t pos, off_t last_pos, int Lseq, int nseq);
void free_sequence(char **seq, int nseq);
int*  find_snps(char **seq, int nseq, int Lseq, int *nsnp, int *nsnp_mult);
char *check_alloc(int nbrelt, int sizelt);
void get_comp(char **seq, int nseq, int Lseq, int *nsite, int *meanLnoN, char *fGC);
void count_rec(char **seq, int ni, int nsnp, int *listSNP, char *nrec, char *snpcov, int verbose, char *locusID);
int test_rec(char **seq, int ni, int snp1, int snp2, int *nhap);
int get_num_non_singleton(int *allele_count);



#include "RNAseqFGT.h"


/*
cc -Wall -o RNAseqFGT RNAseqFGT.c RNAseqFGT_seq_reading.c RNAseqFGT_analysis.c -I RNAseqFGT.h

RNAseqFGT test.fas test.out -v 



*/



int main(int argc, char *argv[])
{
FILE *out_file,  *fasta_file;
char string[MAX_LINE+1];
char ID_locus[MAX_LINE], ID_species[MAX_LINE], ID_individual[MAX_LINE], ID_allele[MAX_LINE];
char ID_previous_contig[MAX_LINE];
int verbose = 0;

int ntot=0;
off_t pos;


if(argc < 3) {
	fprintf(stderr, "%s version 2.0: search for recombination events using the four-gamete test (FGT) on unphased RNAseq data (PopPhyl-like).\nUsage: %s fasta_file out_file [-v] \n", argv[0], argv[0]);
	fprintf(stderr, "\nOption -v (verbose): print on stdout the details on detected SNPs and recombination events. \n");
	exit(1);
	}
	
	
if((fasta_file = fopen(argv[1], "r")) == NULL) {
	fprintf(stderr, "EXIT: cannot open file %s\n", argv[1]);
	exit(1);
	}

if((out_file = fopen(argv[2], "w")) == NULL) {
	fprintf(stderr, "EXIT: cannot create file %s\n", argv[2]);
	exit(1);
	}

if(argc == 4) verbose = 1;

// Check line length
while(fgets(string, MAX_LINE, fasta_file) != NULL) {
	string[MAX_LINE]=0;
	if(strlen(string) >= (MAX_LINE-2)) {
		fprintf(stderr, "EXIT: input file (%s) not in correct FASTA format (line length > %d)\n", argv[1], MAX_LINE);
		exit(1);
		}	
	}
rewind(fasta_file);

// Reads the FASTA file to identify each locus and run its analysis
strcpy(ID_previous_contig, "");
pos = ftello(fasta_file);

while(fgets(string, MAX_LINE, fasta_file) != NULL) {
	
	
	// Get sequence IDs
	if(string[0] == '>') {
		get_IDs(string, ID_locus, ID_species, ID_individual, ID_allele);
		
		// If find a new contig, then run the analysis for this contig
		if(strcmp(ID_previous_contig, ID_locus) != 0) {
			ntot++;
			
			run_analysis(fasta_file, pos, out_file, ID_locus, ID_species, verbose);
			strcpy(ID_previous_contig, ID_locus);
			fseeko(fasta_file, pos, 0);
			}
		
		
		}
	pos = ftello(fasta_file);
	}

fprintf(stderr, "------------\nNormal program end. %d loci analyzed.\n", ntot);


exit(0);
}
/***************************************************************************/


/* Run the analysis for a given contig */

void run_analysis(FILE *fasta_file, off_t pos, FILE *out_file, char *ID_locus, char *ID_species, int verbose)
{
char **seq;    // Aligned sequences
int Lseq;      // Length of aligned sequences
int nsite_noN; // number of sites with at least 2 non-N bases
int meanLnoN;  // mean length of sequences (without N's)
int ni;        // number of individuals 
int nsnp;      // number of bi-allelic SNPs 
int nsnp_mult; // number of SNPs with >2 alleles 
int *listSNP;  // list of positions of bi-allelic SNPs
char fGC[10];  // GC-content
char nrec[10]; // Number of recombination events detected (NA if recombination is not detectable: less than 2 SNPs, or less than 4 identifiable haplotypes)
char snpcov[10];// Average coverage at SNP positions (i.e. proportion of non-N sequences)

off_t last_pos;
static int test = 0;


// Print header
if(test == 0) {
	fprintf(out_file, "LocusID\tSpecies\tN_individuals\tAlnSeqLength\tN_sites_noN\tAvgSeqLg_noN\tN_SNPsBi\tN_SNPs_Mult\tGC_content\tNbRec\tSNP_coverage\n");
	test = 1;

	if(verbose == 1) printf("For each locus, the program performs the four-gamete test on all pairs of SNPs (intervals), to detect recombination events. \
Then the program counts the number of non-overlapping recombinant intervals (to avoid counting twice the same recombination \
event). For each locus, the list of detected SNPs, the total list of recombinant intervals and the subset of non-overlapping \
recombinant intervals are detailed below: \n");
	
	}
	
// Check sequences and count the number of individuals
ni = check_contig_sequences(fasta_file, pos, &last_pos, ID_locus, &Lseq);

// Load sequences
seq = load_sequence(fasta_file, pos, last_pos, Lseq, ni*2);

// compute base composition and sequence length parameters (number of sites with at least 2 non-N bases; mean length of sequences (without N's)
get_comp(seq, ni*2, Lseq, &nsite_noN, &meanLnoN, fGC);

// find SNPs
listSNP = find_snps(seq, ni*2, Lseq, &nsnp, &nsnp_mult);

// Search for recombination events
count_rec(seq, ni, nsnp, listSNP, nrec, snpcov, verbose, ID_locus);

fprintf(out_file, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", ID_locus, ID_species, ni, Lseq, nsite_noN, meanLnoN, nsnp, nsnp_mult, fGC, nrec, snpcov);

//fprintf(out_file, "snp_pos:");
//for(ii=0; ii<nsnp;ii++) fprintf(out_file, " %d", listSNP[ii]);
//fprintf(out_file, "\n");

free_sequence(seq, ni*2);
free(listSNP);


}
/***************************************************************************/





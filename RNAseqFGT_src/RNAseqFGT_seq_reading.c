#include "RNAseqFGT.h"


/**********************************************************/
/* get_IDs: Retrieve all IDs from an ID line              */
/**********************************************************/

void get_IDs(char *string, char *ID_locus, char *ID_species, char *ID_individual, char *ID_allele)
{
int ii, jj;
int n = 0;

jj=0;
for(ii=1;ii<MAX_LINE;ii++) {
	if(string[ii]=='|') {
		ID_locus[jj] = 0;
		n++;
		break;
		}
	ID_locus[jj++] = string[ii];
	}

ii++;
jj=0;
for(;ii<MAX_LINE;ii++) {
	if(string[ii]=='|') {
		ID_species[jj] = 0;
		n++;
		break;
		}
	ID_species[jj++] = string[ii];
	}

ii++;
jj=0;
for(;ii<MAX_LINE;ii++) {
	if(string[ii]=='|') {
		ID_individual[jj] = 0;
		n++;
		break;
		}
	ID_individual[jj++] = string[ii];
	}

ii++;
jj=0;
for(;ii<MAX_LINE;ii++) {
	if(string[ii]=='\n') {
		ID_allele[jj] = 0;
		n++;
		break;
		}
	ID_allele[jj++] = string[ii];
	}


//printf("n=%d %s %s %s %s\n", n, ID_locus, ID_species, ID_individual, ID_allele);
//exit(0);
        
if(n!=4) {
		fprintf(stderr, "EXIT: input file not in correct FASTA format. Expected ID line :\n>ID_locus|ID_species|ID_individual|ID_allele\n");
		exit(1);
	
	}

}
/***************************************************************************/



/************************************************************************/
/* check_contig_sequences: counts the number of individuals for a given */
/* contig + get sequence length + gets the last position in the         */
/* FASTA file + check that each individual has 2 sequences, and that    */
/* all sequences have the same length.                                  */
/************************************************************************/
 


int check_contig_sequences(FILE *fasta_file, off_t first_pos, off_t *last_pos, char *ID_target_contig, int *Lseq)
{
char string[MAX_LINE+1];
char ID_locus[MAX_LINE], ID_species[MAX_LINE], ID_individual[MAX_LINE], ID_allele[MAX_LINE];
int ni = 0, ii;
int test = 2;
char ID_previous_individual[MAX_LINE];
off_t pos;
int lg, lg1;


// Get the last position for that contig in the FASTA file
fseeko(fasta_file, first_pos, 0);
while(fgets(string, MAX_LINE, fasta_file) != NULL) {

	if(string[0] == '>') {
		get_IDs(string, ID_locus, ID_species, ID_individual, ID_allele);
		
		// If not the target contig, then stop
		if(strcmp(ID_target_contig, ID_locus) != 0) {
			break;
			}
		ni++;		
		}
	*last_pos = ftello(fasta_file);
	}




// Check that each individual has 2 alleles
fseeko(fasta_file, first_pos, 0);
strcpy(ID_previous_individual, "");

while(fgets(string, MAX_LINE, fasta_file) != NULL) {
	pos = ftello(fasta_file);
	if(pos > *last_pos) break;
	
	if(string[0] == '>') {
		get_IDs(string, ID_locus, ID_species, ID_individual, ID_allele);
		
		
		// Check that each individual has exactly 2 sequences (alleles)
		if(strcmp(ID_previous_individual, ID_individual) == 0) {
			test++;
			if(test > 2) {
				fprintf(stderr, "EXIT: contig %s: more than two alleles for individual %s\n", ID_locus, ID_individual);
				exit(1);
				}

			}
		else {
			if(test < 2) {
				fprintf(stderr, "EXIT: contig %s: less than two alleles for individual %s\n", ID_locus, ID_previous_individual);
				exit(1);
				}
			test = 1;
			strcpy(ID_previous_individual, ID_individual);

			
			}
		
		}
	}

if(test != 2) {
	fprintf(stderr, "EXIT: contig %s: individual %s does not have 2 alleles\n", ID_locus, ID_previous_individual);
	exit(1);
	}


// Check that all sequences have the same length
fseeko(fasta_file, first_pos, 0);
ii = lg = 0;

while(fgets(string, MAX_LINE, fasta_file) != NULL) {
	pos = ftello(fasta_file);
	if(pos > *last_pos) break;
	
	if(string[0] == '>') {
		ii++;
		get_IDs(string, ID_locus, ID_species, ID_individual, ID_allele);
		
		if(ii > 2 & lg1 != lg) {
			fprintf(stderr, "EXIT: sequence %s|%s|%s|%s (%d bp) does not have the same length as the previous sequence from the same contig (%d bp)\n", ID_locus, ID_species, ID_individual, ID_allele, lg, lg1);
			exit(1);
			}
		lg1 = lg;
		lg = 0;
		}
	else {
		check_sequence(string, ID_locus, ID_species, ID_individual, ID_allele);
		lg = lg + strlen(string);
		}
	}

if(lg1 != lg) {
	fprintf(stderr, "EXIT: sequence %s|%s|%s|%s (%d bp) does not have the same length as the previous sequence from the same contig (%d bp)\n", ID_locus, ID_species, ID_individual, ID_allele, lg, lg1);
	exit(1);
	}


*Lseq = lg;

return(ni/2);
}
/***************************************************************************/


/************************************************************************/
/* check_sequence: change to upper case, check that sequences           */
/* contain only ACGTN bases, and remove newline/blank characters.       */
/************************************************************************/

void check_sequence(char *string, char *ID_locus, char *ID_species, char *ID_individual, char *ID_allele)
{
int ii, jj;
char buf[MAX_LINE];
char c;


jj = 0;
for(ii=0;ii<MAX_LINE; ii++) {
	c = toupper(string[ii]);
	
	if(c == '\n') {buf[jj]=0; break;}
	if(c == ' ' | c == '\t') continue;
	
	// Indels are considered as non-informative sites (N's).
	if(c == '-') c = 'N';
	if(c != 'A' & c != 'C' & c != 'G' & c != 'T' & c != 'N') {
		fprintf(stderr, "EXIT: sequence %s|%s|%s|%s contains non-ACGTN bases (%c)\n", ID_locus, ID_species, ID_individual, ID_allele, string[ii]);
		exit(1);
		}
	
	buf[jj++]=c;
	}

strcpy(string, buf);

}
/***************************************************************************/


/************************************************************************/
/* check_alloc: allocate memory.                                        */
/************************************************************************/
char *check_alloc(int nbrelt, int sizelt)
{
char *retval;
if( (retval=calloc(nbrelt,sizelt)) != NULL ) return(retval);


fprintf(stderr,"\nERROR: not enough memory (check_alloc).\n");


exit(1);
}
/***************************************************************************/


/************************************************************************/
/* load_sequence: read sequences                                        */
/* (NB: to be used only after the check_contig_sequences function.      */
/************************************************************************/

char** load_sequence(FILE *fasta_file, off_t first_pos, off_t last_pos, int Lseq, int nseq)
{
char string[MAX_LINE+1];
char **seq;
int ii;
off_t pos;

seq = (char **) check_alloc(nseq, sizeof(char*));

for(ii = 0; ii< nseq; ii++) {
	seq[ii] = (char *) check_alloc(Lseq+1, sizeof(char));
	strcpy(seq[ii], "");
	}


// Reads all sequences 
fseeko(fasta_file, first_pos, 0);
ii = -1;

while(fgets(string, MAX_LINE, fasta_file) != NULL) {

	pos = ftello(fasta_file);
	if(pos > last_pos) break;
	
	if(string[0] == '>') {
		ii++;
		}
	else {
		check_sequence(string, "na", "na", "na", "na");
		strcat(seq[ii], string);
		}
	}

return(seq);
}
/***************************************************************************/



/************************************************************************/
/* free_sequence: free memory                                           */
/************************************************************************/

void free_sequence(char **seq, int nseq)
{
int ii;

for(ii = 0; ii< nseq; ii++) free(seq[ii]);
free(seq);
}
/***************************************************************************/





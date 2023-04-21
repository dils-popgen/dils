#include "RNAseqFGT.h"



/************************************************************************/
/* find_snps: count the number of SNPs (bi-allelic or with >2 alleles) */
/* and record the location of bi-allelic SNPs                           */
/************************************************************************/

int* find_snps(char **seq, int nseq, int Lseq, int *nsnp, int *nsnp_mult)
{
int ii, jj;
int A, C, G, T, N;
char c;
int sum;
int *listSNP;

listSNP = (int *) check_alloc(Lseq, sizeof(int*));


*nsnp = 0;
*nsnp_mult = 0;


for(ii=0;ii<Lseq;ii++) {
	A = C = G = T = N = 0;
	
	for(jj=0;jj<nseq;jj++) {
		c = seq[jj][ii];
		
		if(c == 'A') A = 1;
		else if(c == 'C') C = 1;
		else if(c == 'G') G = 1;
		else if(c == 'T') T = 1;
		else if(c == 'N') N = 1;
		
		}
	
	sum = A + C + G + T;
	
	if(sum == 2) {
		listSNP[*nsnp] = ii;
		*nsnp = *nsnp + 1;
		}
	else if(sum > 2) *nsnp_mult = *nsnp_mult + 1;
	}

return(listSNP);

}
/***************************************************************************/




/************************************************************************/
/* get_comp: compute the base composition (G+C content), the number of  */
/* aligned sites with at least 2 non-N bases, and the mean length of    */
/* sequences (after excluding all N-sites)                              */
/************************************************************************/

void get_comp(char **seq, int nseq, int Lseq, int *nsite, int *meanLnoN, char *fGC)
{
int ii, jj;
int A, C, G, T, N;
int Acum, Ccum, Gcum, Tcum, Ncum;
char c;
int sum;


*nsite = 0;

Acum = Ccum = Gcum = Tcum = Ncum = 0;

for(ii=0;ii<Lseq;ii++) {
	A = C = G = T = N = 0;
	
	for(jj=0;jj<nseq;jj++) {
		c = seq[jj][ii];
		
		if(c == 'A') A++;
		else if(c == 'C') C++;
		else if(c == 'G') G++;
		else if(c == 'T') T++;
		else if(c == 'N') N++;
		
		}
		
	Acum = Acum + A;
	Ccum = Ccum + C;
	Gcum = Gcum + G;
	Tcum = Tcum + T;
	Ncum = Ncum + N;
	
	sum = A + C + G + T;

	if(sum >= 2) {
		*nsite = *nsite + 1;
		}

	}

sum = Acum + Ccum + Gcum + Tcum;
*meanLnoN = sum / nseq;


if(sum > 0) sprintf(fGC, "%.3f", (float) (Ccum + Gcum)/sum);
else strcpy(fGC, "NA");

}
/***************************************************************************/



/************************************************************************/
/* count_rec: count the number of recombination events detected by the  */
/* four-gamete test. Returns "NA" if recombination is not detectable    */
/* (less than 2 SNPs, or less than four identifiable haplotypes).       */
/* Count also the average coverage at SNP positions (i.e. fraction of   */
/* non-N sequences.                                                     */
/************************************************************************/
void count_rec(char **seq, int ni, int nsnp, int *listSNP, char *nrec, char *snpcov, int verbose, char *locusID)
{
int ii, jj, kk, nn, pos;
int valid_FGT;    // check whether recombination is potentially identifiable
int max_valid_FGT; // check whether recombination is potentially identifiable
int rec, rec_tot, rec_nr;
int n_int, l1, l2;
struct INTERVAL {
	int a;
	int b;
	int na;
	int nb;
	int test;
	} *intervals;
	




// Recode alleles into integer and count number of SNP positions that are non-N
kk = nn = 0;
for(ii = 0; ii < nsnp; ii++) {
	pos = listSNP[ii];
	
	for(jj=0;jj<(2*ni);jj++) {
		kk++;
		if(seq[jj][pos] == 'N') { seq[jj][pos] = 0; nn++;}
		else if(seq[jj][pos] == 'A') seq[jj][pos] = 1;
		else if(seq[jj][pos] == 'C') seq[jj][pos] = 2;
		else if(seq[jj][pos] == 'G') seq[jj][pos] = 3;
		else if(seq[jj][pos] == 'T') seq[jj][pos] = 4;
		else {fprintf(stderr, "EXIT: error (3).\n"); exit(1);}
		}	
	}

if(nsnp==0) strcpy(snpcov, "NA");
else sprintf(snpcov, "%.1f", (float) (kk - nn) / (nsnp * ni * 2));


// If less than 2 SNPs, it is not possible to detect recombination
if(nsnp < 2) {
	if(verbose == 1) {
		printf("---------\n");
		printf("%s\tSNP positions (N=%d):", locusID, nsnp);
		for(ii = 0; ii < nsnp; ii++) {
			printf("\t%d (%d)", listSNP[ii]+1, ii+1);
			}
		printf("\n");
		}

	strcpy(nrec, "NA");
	return;
	}


// Consider all possible pairs of SNPs and run the FGT
n_int = nsnp*(nsnp-1)/2;
intervals = (struct INTERVAL *) check_alloc(n_int, sizeof(struct INTERVAL));

max_valid_FGT = rec_tot = 0;

for(ii = 0; ii < nsnp; ii++) {
	for(jj = ii + 1; jj  < nsnp; jj++) {
		
		rec = test_rec(seq, ni, listSNP[ii], listSNP[jj], &valid_FGT);
		if(valid_FGT > max_valid_FGT) max_valid_FGT = valid_FGT;
		
		// If recombination detected between i and j, then record interval
		if(rec == 1) {
			intervals[rec_tot].a = listSNP[ii];
			intervals[rec_tot].b = listSNP[jj];
			intervals[rec_tot].na = ii + 1;
			intervals[rec_tot].nb = jj + 1;
			intervals[rec_tot].test = 1;
			rec_tot++;
			
			}
		}
	}

// Count the number of non-overlapping intervals
// For each interval, tests if it overlaps with another one,
// and retains only the shortest one.

for(ii = 0; ii < rec_tot; ii++) for(jj = ii+1; jj < rec_tot; jj++) {
	//Ignore intervals that have already been excluded
	if(intervals[ii].test == 0) continue;
	if(intervals[jj].test == 0) continue;
	
	//Ignore non-overlapping intervals
	if(intervals[ii].b <= intervals[jj].a) continue;
	
	// Overlapping intervals: retain the shortest one
	l1 = intervals[ii].b - intervals[ii].a;
	l2 = intervals[jj].b - intervals[jj].a;
	if(l1 > l2) intervals[ii].test = 0;
	else intervals[jj].test = 0;
	
	}

rec_nr = 0;
for(ii = 0; ii < rec_tot; ii++) if(intervals[ii].test == 1) rec_nr++;

// This case should not occur
if(rec_tot > 0 & rec_nr == 0) {
	fprintf(stderr,"EXIT: error (2) in the four-gamete test.\n");
	exit(1);
	}

if(verbose == 1) {
	printf("---------\n");
	printf("%s\tSNP positions (N=%d):", locusID, nsnp);
	for(ii = 0; ii < nsnp; ii++) {
		printf("\t%d (%d)", listSNP[ii]+1, ii+1);
		}
	printf("\n");
	}
if(verbose == 1 & rec_nr > 0) {

	printf("%s\tAll detected recombinant intervals:", locusID);
	for(ii = 0; ii < rec_tot; ii++) {
		printf("\t%d..%d (%d..%d)", intervals[ii].a, intervals[ii].b, intervals[ii].na, intervals[ii].nb);
		}
	printf("\n");
	
	printf("%s\tNon-overlapping recombinant intervals:", locusID);
	for(ii = 0; ii < rec_tot; ii++) if(intervals[ii].test == 1) {
		printf("\t%d..%d (%d..%d)", intervals[ii].a, intervals[ii].b, intervals[ii].na, intervals[ii].nb);
		}
	printf("\n");
	}


free(intervals);

// To be able to detect recombination, at least one pair of non-singleton SNPs with four identifiable gametes is required.
// If this criteria is not met, then return NA
if(max_valid_FGT == 0) {
	strcpy(nrec, "NA");
	return;
	}
	
sprintf(nrec, "%d", rec_nr);
}
/***************************************************************************/





/************************************************************************/
/* test_rec: identify haplotypes and performs the four-gamete test for  */
/* a pair of SNPs.                                                      */
/* returns the number of recombination events detected                  */
/* valid_FGT is set to 1 if recombination is potentially detectable, O  */
/* otherwise */
/************************************************************************/
int test_rec(char **seq, int ni, int snp1, int snp2, int *valid_FGT)
{
int haplotypes[5][5]; // count matrix of all possibles haplotypes (A,C,G,T,N) x (A,C,G,T,N)
// For each SNP, count number of occurence of each allele among identifiable haplotypes:
int allele_count_snp1[5], allele_count_snp2[5]; 

int ii, jj, kk, i1, i2, n1, n2;
int a1, a2; // the 2 alleles at the 1st SNP
int b1, b2; // the 2 alleles at the 2nd SNP

int ngam;    // number of identifiable gametes

// Initiate count matrix of all possibles haplotypes
for(ii=0;ii<5;ii++) for(jj=0;jj<5;jj++) haplotypes[ii][jj] = 0;

// Initiate allele counts
for(ii=0;ii<5;ii++) allele_count_snp1[ii] = allele_count_snp2[ii] = 0;

ngam = 0;

// For each individual, identify haplotypes (when possible)
for(ii = 0; ii < ni; ii++) {
	
	i1 = 2 * ii;
	i2 = i1 + 1;
	
	a1 = seq[i1][snp1];
	a2 = seq[i2][snp1];
	
	b1 = seq[i1][snp2];
	b2 = seq[i2][snp2];
	
	
	// If both SNPs are heterozygotes then it is not possible to infer haplotypes
	if(a1 != a2 & b1 != b2) continue;
		
	else {
		// For each SNP, count number of occurence of each allele among identifiable haplotypes
		allele_count_snp1[a1] = allele_count_snp1[a1] + 1;
		allele_count_snp1[a2] = allele_count_snp1[a2] + 1;
		allele_count_snp2[b1] = allele_count_snp2[b1] + 1;	
		allele_count_snp2[b2] = allele_count_snp2[b2] + 1;	

	 	// If both SNPs are homozygote then we genotype one gamete (haplotype)
		if(a1 == a2 & b1 == b2) {
			haplotypes[a1][b1] =  1;
			ngam = ngam + 1;
			}

		// Otherwise (1 homozygote + 1 heterozygote) then we genotype two gametes (haplotypes)
		else {
			haplotypes[a1][b1] =  1;
			haplotypes[a2][b2] =  1;
			ngam = ngam + 2;
			}
		}
	}

// Count the number of (non-N) different haplotypes

kk = 0;
for(ii=1;ii<5;ii++) for(jj=1;jj<5;jj++) kk = kk + haplotypes[ii][jj];

// Count the number of (non-N) non-singleton alleles for each SNP
n1 = get_num_non_singleton(allele_count_snp1);
n2 = get_num_non_singleton(allele_count_snp2);


// by default valid_FGT is set to 1
*valid_FGT = 1;

// If one of the alleles is a singleton, then it is impossible to detect recombination
if(n1 < 2 || n2 < 2) *valid_FGT = 0;

// If less than 4 identifiable gametes, then it is impossible to detect recombination
if(ngam < 4) *valid_FGT = 0;



// If 4 different haplotypes => recombination is detected
if(kk == 4) return(1);
// Otherwise, no recombination is detected
if(kk < 4) return(0);
// Other cases should not occur
fprintf(stderr,"EXIT: error in the four-gamete test.\n");
exit(1);

}
/***************************************************************************/


// Count the number of (non-N) non-singleton allele 

int get_num_non_singleton(int *allele_count)
{
int ii;
int nb_non_singleton=0;

for(ii=1;ii<5;ii++) if(allele_count[ii] > 1) nb_non_singleton++;

return(nb_non_singleton);

}
/***************************************************************************/













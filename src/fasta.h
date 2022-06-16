
typedef struct
{
	char *name;
	char *name2;
	char *path;
	char *SCFname;
	int  length;
	char *data;
	char *qual;
	int  finished;
} fasta;

#define B64_long long int
fasta *decodeFastq (char *fname, int *nContigs, long *tB, char* pdata, long Size_pdata,fasta *segg);
void fastaLC (fasta *seg, int nSeg);
void fastaUC (fasta *seg, int nSeg);


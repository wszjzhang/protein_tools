#ifndef TMSCORE_H
#define TMSCORE_H

#define  NMAX 1500

static char aa[3*22+1] = "BCKGLYALASERCYSVALTHRILEPROMETASPASNLEULYSGLUGLNARGHISPHETYRTRPCYX";
static char slc[1*22+1] = "XGASCVTIPMDNLKEQRHFYWC";	


typedef struct  {	
//for both tm-score and tm-align of one alignment
	char    pcProtID[20];
	int     nSeqLen;
	char    pcAliQ[NMAX];
	char    pcAliS[NMAX];
	
	double  lfRmsd;
	double  lfTMscore;
	double  lfGDT;
	double  lfGDTHA;	
	
	int     nHitLen;
	int     nAliLen;	
	int     nIdenNum;
	int     nGapNum;
	int     nPosNum;	
	int     nStartQ;
	int     nStartS;
	double  lfPsibScore;	
	double  zscore;
}TSPScore;



typedef struct
{

	double xa[NMAX];
	double ya[NMAX];
	double za[NMAX];


	double xt[NMAX];
	double yt[NMAX];
	double zt[NMAX];


	double xb[NMAX];
	double yb[NMAX];
	double zb[NMAX];
}StruInfo;

/*typedef struct
{

}NresInfo;*/

typedef struct
{
	char seqA[NMAX];
	char seqB[NMAX];
	char seq1A[NMAX];
	char seq1B[NMAX];

	int nresA[NMAX];
	int nresB[NMAX];
	int nseqA;
	int nseqB;
}SeqInfo;


typedef struct
{
	double d;
	double d0;
	double d0_fix;

	double d0_search;
	double d_output;
}DPara;


typedef struct
{
	int n_ali;
	int iA[NMAX];
	int iB[NMAX];

	int i_ali[NMAX];
	int n_cut; //[1,n_ali],align residues for the score

	int ka0;
	int k_ali[NMAX];
	int k_ali0[NMAX];
}AlignInfo;


typedef struct
{
	double score;
	double score_max;       //TM-score

	double score_fix;
	double score_fix_max;

	double score_maxsub;	//MaxSub-score
	double score_maxsub_max;

	double score10;         //TM-score10
	double score10_max;	
}ScoresInfo;


typedef struct
{
	int n_GDT05;
	int n_GDT1;
	int n_GDT2;
	int n_GDT4;
	int n_GDT8;

	int n_GDT05_max;            //!number of residues<0.5
	int n_GDT1_max;             //!number of residues<1
	int n_GDT2_max;             //!number of residues<2
	int n_GDT4_max;             //!number of residues<4
	int n_GDT8_max;             //!number of residues<8
}GdtInfo;


typedef struct
{
	int m_out;
	int m_fix;
	char outname[100];
	char pdb[2][100];	

	int  m_outResult;
	char resultFileName[100];
	//char testFileName[100];
}InputPara;


typedef struct
{
	int n_it;
	int n_init;	
	int L_ini[100];
}IterateInfo;


typedef struct
{	
	double rmsd;
	double armsd;
	double rmsd_ali;
}RmsdInfo;

int TMscore(StruInfo &stru, SeqInfo &seqs, AlignInfo &aligns, ScoresInfo &scores, GdtInfo &gdt, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, double plfSuperMatrix[3][4]);
int programHelpInfo(int argcNum, char **argvPara, InputPara &inputPara);
int getInputPara(int argcNum, char **argvPara, InputPara &inputPara, DPara &dPara);
int readFirstCaFile(char *pdb, StruInfo &stru,  SeqInfo &seqs);
int readSecondCaFile(char *pdb, StruInfo &stru,  SeqInfo &seqs);
int pickupAlignedRes( AlignInfo &aligns, SeqInfo &seqs, StruInfo &stru);
int generateParameter( AlignInfo &aligns,SeqInfo &seqs, InputPara &inputPara, DPara &dPara, IterateInfo &iterate);
int findMaxScore(StruInfo &stru, SeqInfo &seqs, AlignInfo &aligns, DPara &dPara, ScoresInfo &scores, GdtInfo &gdt, IterateInfo &iterate, RmsdInfo &rmsds);
int	outputResults(StruInfo &stru, SeqInfo &seqs, AlignInfo &aligns, ScoresInfo &scores, GdtInfo &gdt, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, double plfSuperMatrix[3][4]);
int outputRotatedChain(StruInfo &stru,  AlignInfo &aligns, ScoresInfo &scores, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, SeqInfo &seqs);	
int outputAlignedSequences(StruInfo &stru,  AlignInfo &aligns, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, SeqInfo &seqs);
int u3b(double w[NMAX], double x[3][NMAX], double y[3][NMAX], int n, int mode, double &rms, double u[3][3], double t[3], int &ier);
int score_fun(StruInfo &stru, SeqInfo &seqs, AlignInfo &aligns, DPara &dPara, ScoresInfo &scores, GdtInfo &gdt);

void transSeq2Seq3(char *ss1, char *ss3, int nSeqLen);
#endif

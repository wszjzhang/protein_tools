/*********************************************************************************
*                 Copyright (c) 2007 Digital Biology Laboratory
*                           All Rights Reserved
*
*  FILE:     calGdtRmsdScore.cpp
*  PURPOSE:
*            This program is to calculate the GDT_TS and Rmsd score of two alignned 
*			 structures
*
*  This is unpublished proprietary source code.
*  Date:
*  PI:       Dong Xu
*  AUTHOR:   Jingfen Zhang
*            
*
*  NOTE:     
*
*  CHANGE HISTORY:
*            First version: 2008-04
*
**********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "TMSM.h"

#define MAXPATHLEN 500
#define MAXSEQLEN 1500

void calTMscore_PD(char *pcPdbFName, char *pcDBFName, double &lfRmsd, double &lfGdt, double &lfGdtHA, double &lfTmscore, double plfSuperMatrix[3][4]);
void calTMscore_PP(char *pcPdbFName1, char *pcPdbFName2, double &lfRmsd, double &lfGdt, double &lfGdtHA, double &lfTmscore, double plfSuperMatrix[3][4]);
void getDBCord(char *pcDBFName, int &nSeqLen);
void getDBCord(char *pcDBFName, double *plfDBCord, char *pcDBSeq, int &nSeqLen);
void getPdbCord(char *pcPdbFName, double *plfPdbCord, char *pcPdbSeq);
void calSeqAlign_TMScore(double *plfCaCord1, double *plfCaCord2, char *pcAli, char *pcProtID, double &lfRmsd, double &lfTmscore, double &lfGdtHA, double &lfGdt, double plfSuperMatrix[3][4]);
void AAP_setStruSeqInfo(StruInfo &stru, SeqInfo &seqs, double *plfCaCord1, double *plfCaCord2, char *pcAli);
void transSeq3Seq(char *ss1, char *ss3, int nSeqLen);
void chompx(char *str);
int  AAP_pickupAlignedRes(AlignInfo &aligns, SeqInfo &seqs, StruInfo &stru);

/***********************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
***********************************************************************/
int 
main(int argc, char **argv)
{
	char pcPdbFName1[MAXPATHLEN], pcPdbFName2[MAXPATHLEN];
	double lfRmsd, lfGdt, lfGdtHA, lfTmscore, plfSuperMatrix[3][4];

	if(argc == 1)
	{
		fprintf(stderr, "calGdtRmsdScore PDBFName1 PDBFName2(or DBFName)\n");
		exit(0);
	}
	
	strcpy(pcPdbFName1, argv[1]);
	strcpy(pcPdbFName2, argv[2]);
	//calTMscore_PP(pcPdbFName1, pcPdbFName2, lfRmsd, lfGdt, lfGdtHA, lfTmscore, plfSuperMatrix);
	//fprintf(stderr,"%-40s %-40s %-8.3f %-8.3f %-8.3f %-8.3f\n", pcPdbFName1, pcPdbFName2, lfRmsd, lfGdt, lfGdtHA, lfTmscore, plfSuperMatrix);

	calTMscore_PD(pcPdbFName1, pcPdbFName2, lfRmsd, lfGdt, lfGdtHA, lfTmscore, plfSuperMatrix);
	fprintf(stdout,"rmsd gdt gdtHA TMscore %.3f %.3f %.3f %.3f\n", lfRmsd, lfGdt, lfGdtHA, lfTmscore);
	
	return 0;
}

/***********************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
***********************************************************************/
void 
calTMscore_PD(char *pcPdbFName, char *pcDBFName, double &lfRmsd, double &lfGdt, double &lfGdtHA, double &lfTmscore, double plfSuperMatrix[3][4])
{
	char   pcPdbSeq[MAXSEQLEN], pcDBSeq[MAXSEQLEN], pcAli[MAXSEQLEN];
	double plfPdbCord[MAXSEQLEN*3], plfDBCord[MAXSEQLEN*3];
	int nAANo, nSeqLen;

	getDBCord(pcDBFName, plfDBCord, pcDBSeq, nSeqLen);
	getPdbCord(pcPdbFName, plfPdbCord, pcPdbSeq);

	//fprintf(stderr,"seqLen = %d %s\n", nSeqLen, pcDBSeq);
	//fprintf(stderr,"%s\n", pcPdbSeq);
	memset(pcAli, '\0', MAXSEQLEN);
	for(nAANo = 0; nAANo < nSeqLen; nAANo++)
	{
		if(pcDBSeq[nAANo] == '-' || pcPdbSeq[nAANo] == '-')
			pcAli[nAANo] = '-';
		else
			pcAli[nAANo] = pcDBSeq[nAANo];
	}
	
	calSeqAlign_TMScore(plfPdbCord, plfDBCord, pcAli, "test", lfRmsd, lfTmscore, lfGdtHA, lfGdt, plfSuperMatrix);
	//fprintf(stderr,"%-8.3f %-8.3f %-8.3f %-8.3f\n", lfTmscore,  lfGdt, lfGdtHA, lfRmsd);
	return;
}


/***********************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
***********************************************************************/
void 
calTMscore_PP(char *pcPdbFName1, char *pcPdbFName2, double &lfRmsd, double &lfGdt, double &lfGdtHA, double &lfTmscore, double plfSuperMatrix[3][4])
{
	char   pcPdbSeq1[MAXSEQLEN], pcPdbSeq2[MAXSEQLEN], pcAli[MAXSEQLEN];
	double plfPdbCord1[MAXSEQLEN*3], plfPdbCord2[MAXSEQLEN*3];
	int    nSeqLen, nAANo;

	getPdbCord(pcPdbFName1, plfPdbCord1, pcPdbSeq1);
	getPdbCord(pcPdbFName2, plfPdbCord2, pcPdbSeq2);

	nSeqLen = strlen(pcPdbSeq1);
	if (nSeqLen > strlen(pcPdbSeq2))
		nSeqLen = strlen(pcPdbSeq2);
	pcPdbSeq1[nSeqLen] = '\0';
	pcPdbSeq2[nSeqLen] = '\0';
	

	memset(pcAli, '\0', MAXSEQLEN);
	for(nAANo = 0; nAANo < nSeqLen; nAANo++)
	{
		if(pcPdbSeq1[nAANo] == '-' || pcPdbSeq2[nAANo] == '-')
			pcAli[nAANo] = '-';
		else
			pcAli[nAANo] = pcPdbSeq1[nAANo];
	}
	
	calSeqAlign_TMScore(plfPdbCord1, plfPdbCord2, pcAli, "test", lfRmsd, lfTmscore, lfGdtHA, lfGdt, plfSuperMatrix);
	return;
}

/***********************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
***********************************************************************/
void
getPdbCord(char *pcPdbFName, double *plfPdbCord, char *pcPdbSeq)
{
	FILE *fp;
	char pcLine[300], pcAANo[7], pcAA[4], pcCord[9];
	int  nCordNo, nAANo;
	memset(plfPdbCord, '\0', MAXSEQLEN*3*sizeof(double));
	memset(pcPdbSeq, '-', (MAXSEQLEN-1)*sizeof(char));
	pcPdbSeq[MAXSEQLEN-1] = '\0';
	fp = fopen(pcPdbFName, "r");
	while(!feof(fp))
	{
		memset(pcLine, '\0', 300);
		fgets(pcLine, 300, fp);
		if(strstr(pcLine, "TER"))
			break;
		
		if(strstr(pcLine, "ATOM") && strstr(pcLine, "CA  "))
		{
			memset(pcAANo, '\0', 7);
			strncpy(pcAANo, pcLine+22, 6);
			nAANo = atoi(pcAANo) - 1;
			
			for(nCordNo=0; nCordNo<3; nCordNo++)
			{
				memset(pcCord, '\0', 9);
				strncpy(pcCord, pcLine+30+nCordNo*8, 8);
				sscanf(pcCord, "%lf", plfPdbCord+nAANo*3+nCordNo);				
			}	
			
			memset(pcAA, '\0', 4);
			strncpy(pcAA, pcLine+17, 3);			
			
			if(strcmp(pcAA, "ALA") == 0) {pcPdbSeq[nAANo]='A'; continue;} 
			if(strcmp(pcAA, "ARG") == 0) {pcPdbSeq[nAANo]='R'; continue;} 
			if(strcmp(pcAA, "ASN") == 0) {pcPdbSeq[nAANo]='N'; continue;} 
			if(strcmp(pcAA, "ASP") == 0) {pcPdbSeq[nAANo]='D'; continue;} 
			if(strcmp(pcAA, "ASX") == 0) {pcPdbSeq[nAANo]='B'; continue;} 
			if(strcmp(pcAA, "CYS") == 0) {pcPdbSeq[nAANo]='C'; continue;} 
			if(strcmp(pcAA, "GLN") == 0) {pcPdbSeq[nAANo]='Q'; continue;} 
			if(strcmp(pcAA, "GLU") == 0) {pcPdbSeq[nAANo]='E'; continue;} 
			if(strcmp(pcAA, "GLX") == 0) {pcPdbSeq[nAANo]='Z'; continue;} 
			if(strcmp(pcAA, "GLY") == 0) {pcPdbSeq[nAANo]='G'; continue;} 
			if(strcmp(pcAA, "HIS") == 0) {pcPdbSeq[nAANo]='H'; continue;} 
			if(strcmp(pcAA, "ILE") == 0) {pcPdbSeq[nAANo]='I'; continue;} 
			if(strcmp(pcAA, "LEU") == 0) {pcPdbSeq[nAANo]='L'; continue;} 
			if(strcmp(pcAA, "LYS") == 0) {pcPdbSeq[nAANo]='K'; continue;} 
			if(strcmp(pcAA, "MET") == 0) {pcPdbSeq[nAANo]='M'; continue;} 
			if(strcmp(pcAA, "PHE") == 0) {pcPdbSeq[nAANo]='F'; continue;} 
			if(strcmp(pcAA, "PRO") == 0) {pcPdbSeq[nAANo]='P'; continue;} 
			if(strcmp(pcAA, "SER") == 0) {pcPdbSeq[nAANo]='S'; continue;} 
			if(strcmp(pcAA, "THR") == 0) {pcPdbSeq[nAANo]='T'; continue;} 
			if(strcmp(pcAA, "TRP") == 0) {pcPdbSeq[nAANo]='W'; continue;} 
			if(strcmp(pcAA, "TYR") == 0) {pcPdbSeq[nAANo]='Y'; continue;} 
			if(strcmp(pcAA, "VAL") == 0) {pcPdbSeq[nAANo]='V'; continue;} 				
		}
		
	}
	fclose(fp);

	for(nAANo = MAXSEQLEN-2; nAANo>0; nAANo--)
	{
		if(pcPdbSeq[nAANo] != '-')
		{
			pcPdbSeq[nAANo+1] = '\0';
			break;
		}
	}

	return;
}

/***********************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
***********************************************************************/
void
getDBCord(char *pcDBFName, double *plfDBCord, char *pcDBSeq, int &nSeqLen)
{
	FILE *fp;
	char pcLine[10*MAXSEQLEN];
	int  nAANo, nLineNo;
	memset(plfDBCord, '\0', MAXSEQLEN*3*sizeof(double));
	memset(pcDBSeq, '\0', MAXSEQLEN*sizeof(char));

	if( (fp = fopen(pcDBFName, "r"))==NULL)
	{
		fprintf(stderr,"%s open error, quit\n", pcDBFName);
		exit(1);
	}
	while(!feof(fp))
	{
		memset(pcLine, '\0', 10*MAXSEQLEN);
		fgets(pcLine, 10*MAXSEQLEN, fp);

		if(strstr(pcLine, ">Real Sequence Info:"))
		{
			memset(pcDBSeq, '\0', MAXSEQLEN);
			fgets(pcDBSeq, MAXSEQLEN, fp);
			chompx(pcDBSeq);
			nSeqLen = strlen(pcDBSeq);		
		}


		if(strstr(pcLine, ">Ca XYZ:"))
		{			
			for(nLineNo = 0; nLineNo<3; nLineNo++)
			{
				memset(pcLine, '\0', 10*MAXSEQLEN);			
				fgets(pcLine, 10*MAXSEQLEN, fp);
				for(nAANo = 0; nAANo<nSeqLen; nAANo++)
					sscanf(pcLine+10*nAANo,"%lf", plfDBCord+3*nAANo+nLineNo); 				    	
			}
			break; 
		}
	}
	fclose(fp);
}

/***********************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
***********************************************************************/
void
getDBCord(char *pcDBFName, int &nSeqLen)
{
	FILE *fp;
	char pcLine[MAXSEQLEN];
	int  nAANo, nLineNo;

	fp = fopen(pcDBFName, "r");
	while(!feof(fp))
	{
		memset(pcLine, '\0', MAXSEQLEN);
		fgets(pcLine, MAXSEQLEN, fp);
		if(strstr(pcLine, ">Real Sequence Info:"))
		{
			memset(pcLine, '\0', MAXSEQLEN);
			fgets(pcLine, MAXSEQLEN, fp);
			chompx(pcLine);
			nSeqLen = strlen(pcLine);		
			break;
		}
	}
	fclose(fp);
}


/***********************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
***********************************************************************/
void 
calSeqAlign_TMScore(double *plfCaCord1, double *plfCaCord2, char *pcAli, char *pcProtID, double &lfRmsd, double &lfTmscore, double &lfGdtHA, double &lfGdt, double plfSuperMatrix[3][4])
{		
	StruInfo    stru;
	SeqInfo     seqs;
	RmsdInfo    rmsds;
	ScoresInfo  scores;
	GdtInfo     gdt;
	AlignInfo   aligns;	
	InputPara   inputPara;
	DPara       dPara;
	
	memset(&stru, '\0', sizeof(StruInfo));	
	memset(&seqs, '\0', sizeof(SeqInfo));
	memset(&scores, '\0', sizeof(ScoresInfo));
	memset(&gdt, '\0', sizeof(GdtInfo));		
	memset(&rmsds, '\0', sizeof(RmsdInfo));	
	memset(&aligns, '\0', sizeof(AlignInfo));
	memset(&inputPara, '\0', sizeof(InputPara));	
	memset(&dPara, '\0', sizeof(DPara));
	
	
	AAP_setStruSeqInfo(stru, seqs, plfCaCord1, plfCaCord2, pcAli);			
	AAP_pickupAlignedRes(aligns, seqs, stru);				
	inputPara.m_out = -1;
	inputPara.m_fix = -1;	
	inputPara.m_outResult = -1;	
		
	strcpy(inputPara.pdb[0], pcProtID);		
	strcpy(inputPara.pdb[1], pcProtID);
	TMscore(stru, seqs, aligns, scores, gdt, inputPara, dPara, rmsds, plfSuperMatrix);	
	
	lfRmsd = rmsds.rmsd_ali;
	lfTmscore = scores.score_max;
	lfGdtHA = 100*(gdt.n_GDT05_max + gdt.n_GDT1_max + gdt.n_GDT2_max + gdt.n_GDT4_max)/(double)(4*seqs.nseqB);
	lfGdt = 100*(gdt.n_GDT1_max + gdt.n_GDT2_max + gdt.n_GDT4_max + gdt.n_GDT8_max)/(double)(4*seqs.nseqB);
	
	return;
}
	

/*********************************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
**********************************************************************************/
void
AAP_setStruSeqInfo(StruInfo &stru, SeqInfo &seqs, double *plfCaCord1, double *plfCaCord2, char *pcAli)
{
	int nAANo, nStruAANum;

	nStruAANum = 0;
	for(nAANo=0; nAANo<strlen(pcAli); nAANo++)
	{
		if(pcAli[nAANo] != '-')
		{			
			*(stru.xa+nStruAANum) = plfCaCord1[3*nAANo+0];
			*(stru.ya+nStruAANum) = plfCaCord1[3*nAANo+1];
			*(stru.za+nStruAANum) = plfCaCord1[3*nAANo+2];
			
			seqs.nresA[nStruAANum] = nStruAANum+1;			
			seqs.seq1A[nStruAANum] = pcAli[nAANo]; 
			
			*(stru.xb+nStruAANum) = plfCaCord2[3*nAANo+0];
			*(stru.yb+nStruAANum) = plfCaCord2[3*nAANo+1];
			*(stru.zb+nStruAANum) = plfCaCord2[3*nAANo+2];
			
			seqs.seq1B[nStruAANum] = pcAli[nAANo];
			seqs.nresB[nStruAANum] = nStruAANum+1;
			
			nStruAANum ++;			
		}
	}

	transSeq2Seq3(seqs.seq1A, seqs.seqA, nStruAANum);	
	seqs.nseqA = nStruAANum;
	transSeq2Seq3(seqs.seq1B, seqs.seqB, nStruAANum);	
	seqs.nseqB = nStruAANum;
	return;
}

/*********************************************************************************
*          NAME:  P2D_genPred2D()
*   DESCRIPTION:  The scheme to evoke the prediction routines.
*         INPUT:  None
*        OUTPUT:  None
*        RETURN:  None
**********************************************************************************/
void 
transSeq3Seq(char *ss1, char *ss3, int nSeqLen)
{
	int nAANo;
	strcpy(ss3, "");
	for(nAANo=0; nAANo<nSeqLen; nAANo++)
	{
		switch(ss1[nAANo])
		{
		case 'G': strcat(ss3,"GLY"); break;
		case 'A': strcat(ss3,"ALA"); break;
		case 'S': strcat(ss3,"SER"); break;
		case 'C': strcat(ss3,"CYS"); break;
		case 'V': strcat(ss3,"VAL"); break;
		case 'T': strcat(ss3,"THR"); break;
		case 'I': strcat(ss3,"ILE"); break;
		case 'P': strcat(ss3,"PRO"); break;
		case 'M': strcat(ss3,"MET"); break;
		case 'D': strcat(ss3,"ASP"); break;
		case 'N': strcat(ss3,"ASN"); break;
		case 'L': strcat(ss3,"LEU"); break;
		case 'K': strcat(ss3,"LYS"); break;
		case 'E': strcat(ss3,"GLU"); break;
		case 'Q': strcat(ss3,"GLN"); break;			
		case 'R': strcat(ss3,"ARG"); break;			
		case 'H': strcat(ss3,"HIS"); break;
		case 'F': strcat(ss3,"PHE"); break;
		case 'Y': strcat(ss3,"TYR"); break;
		case 'W': strcat(ss3,"TRP"); break;			
		case 'X': strcat(ss3,"BCK"); break;			
		case 'Z': strcat(ss3,"GLX"); break;			
		}
	}
	return;
}

/*********************************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
**********************************************************************************/
int
AAP_pickupAlignedRes(AlignInfo &aligns, SeqInfo &seqs, StruInfo &stru)
{
	int i, j, k;	
	k = 0;

	char pcAliQuery[1500], pcAliHit[1500];
	memset(pcAliQuery, '\0', 1500);
	memset(pcAliHit, '\0', 1500);

	for (i = 1; i <= seqs.nseqA; ++i)
	{
		for (j = 1; j <= seqs.nseqB; ++j)
		{
			if (seqs.nresA[i-1] == seqs.nresB[j-1])
			{
				++k;
				aligns.iA[k-1] = i;
				aligns.iB[k-1] = j;

				pcAliHit[k-1] = seqs.seq1A[i-1];				
				pcAliQuery[k-1] = seqs.seq1B[j-1];

				break;
			}
		}
	}	

	aligns.n_ali = k; /* number of aligned residues */
	if (aligns.n_ali < 1)
	{
		printf("There is no common residues in the input structures.\n");
		return -1;
	}
	return 1;
}


/*********************************************************************************
*          NAME:
*   DESCRIPTION:
*         INPUT:
*        OUTPUT:
*        RETURN:
**********************************************************************************/
void
chompx(char *str)
{
	int i;
	for(i=0;(unsigned int)i<strlen(str);i++)
	{
		if(str[i] == 0x0D || str[i] == 0x0A)
			str[i] = 0x0;
	}
	while(1)
	{
		if(str[strlen(str)-1] == 0x20)
			str[strlen(str)-1] = 0x0;
		else
			break;
	}
}

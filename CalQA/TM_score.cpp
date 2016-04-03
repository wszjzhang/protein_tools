/*************************************************************************
*     This program is to compare two protein structures and identify the
*     best superposition that has the maximum TM-score. The program can
*     be freely copied or modified.
*     For comments, please email to: yzhang@ku.edu
*
*     Reference:
*     Yang Zhang, Jeffrey Skolnick, Proteins, 2004 57:702-10.
*
*     Note:
*		The source code have been transfered from Fortune to Standard C by
*       Jingfen Zhang 
******************* Updating history ************************************
*     2005/10/19: the program was reformed so that the score values
*                 are not dependent on the specific compilers.
*     2006/06/20: select 'A' if there is altLoc when reading PDB file
*     2007/02/05: fix a bug with length<15 in TMscore_32
*     2007/02/27: rotation matrix from Chain-1 to Chain-2 is added
*     2007/12/06: GDT-HA score is added, fixed a bug for reading PDB*
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "TMSM.h"


//int main(int argc, char **argv)
int 
TMscore(StruInfo &stru, SeqInfo &seqs, AlignInfo &aligns, ScoresInfo &scores, GdtInfo &gdt, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, double plfSuperMatrix[3][4])
{

	IterateInfo iterate;
	memset(&iterate, '\0', sizeof(IterateInfo));

	/*if(programHelpInfo(argc, argv, inputPara) == 0)
		return -1;

	getInputPara(argc, argv, inputPara);

	if (readFirstCaFile(inputpara.pdb[0], stru, nres, seqs) == -1 )
		return -1;

	if (readSecondCaFile(inputpara.pdb[1], stru, nres, seqs) == -1)
		return -1;	

	if (pickupAlignedRes(aligns, seqs, stru) == -1)
		return -1;
*/	

	generateParameter(aligns, seqs, inputPara, dPara, iterate);
	
	findMaxScore(stru, seqs, aligns, dPara, scores, gdt, iterate, rmsds);

	outputResults(stru, seqs, aligns, scores, gdt, inputPara, dPara, rmsds, plfSuperMatrix);

	if (inputPara.m_out == 1)
		outputRotatedChain(stru, aligns, scores, inputPara, dPara, rmsds, seqs);
	
	outputAlignedSequences(stru, aligns, inputPara, dPara, rmsds, seqs);

	return 0;
}


int programHelpInfo(int argcNum, char **argvPara, InputPara &inputPara)
{
	
	if(argcNum == 1 || strcmp(argvPara[1], "?") == 0 || strcmp(argvPara[1], "-h") == 0)
	{
		fprintf(stderr,"Brief instruction for running TM-score program:\n");
		fprintf(stderr,"1. Run TM-score to compare ''model'' and ''native'':'\n");
		fprintf(stderr,"  TMscore model native\n");
		fprintf(stderr,"2. Run TM-score with an assigned d0, e.g. 5 Angstroms:'\n");
		fprintf(stderr,"  TMscore model native -d 5'\n");
		fprintf(stderr,"3. Run TM-score with superposition output, e.g. TM.sup:'\n");
		fprintf(stderr,"  TMscore model native -o TM.sup'\n");
		fprintf(stderr,"4, To view the superimposed structures by rasmol:'\n");
		fprintf(stderr,"   rasmol -script TM.sup'\n");		
		return 0;
	}
	else
		return 1;
	
}
	
int getInputPara(int argcNum, char **argvPara, InputPara &inputPara, DPara &dPara)
{
	int  nargNo = 0, i=0;
	
// options
	inputPara.m_out=-1;  //decided output
	inputPara.m_fix=-1;  //fixed length-scale only for output
	i=0;
	while(i < argcNum-1)
	{
		i=i+1;
		if(strcmp(argvPara[i], "-o")==0)
		{
			inputPara.m_out=1;
			i = i + 1;
			strcpy(inputPara.outname, argvPara[i]);
		}
		else if(strcmp(argvPara[i], "-d") == 0)
		{
			inputPara.m_fix = 1;
			i = i + 1;
			dPara.d0_fix = atof(argvPara[i]);
		}
		else if(strcmp(argvPara[i], "-R") == 0)
		{
			inputPara.m_outResult = 1;
			i = i + 1;
			strcpy(inputPara.resultFileName, argvPara[i]);
		}
		else
		{
			strcpy(inputPara.pdb[0], argvPara[i]);
			i = i+1;
			strcpy(inputPara.pdb[1], argvPara[i]);			
		}
	}
	return 0;
}
	



int generateParameter( AlignInfo &aligns, SeqInfo &seqs, InputPara &inputPara, DPara &dPara, IterateInfo &iterate)
{
	int i, n_init_max, L_ini_min;

	if(seqs.nseqB > 15)
		dPara.d0 = (double)(pow((double)(seqs.nseqB - 15), (1.0/3.0))*1.24f - 1.8f);
	else
		dPara.d0 = 0.5;
	
	if(dPara.d0 < 0.5)
		dPara.d0 = 0.5;
	
	if(inputPara.m_fix == 1)
		dPara.d0 = dPara.d0_fix;
	
	//d0_search
	dPara.d0_search = dPara.d0;
	
	if(dPara.d0_search > 8)
		dPara.d0_search = 8;
	
	if(dPara.d0_search < 4.5)
		dPara.d0_search = 4.5;
	
//iterative parameters 	
	iterate.n_it = 20;                   //!maximum number of iterations

	dPara.d_output = 4;                //!for output alignment
	if(inputPara.m_fix == 1)
		dPara.d_output = dPara.d0_fix;
	
	n_init_max = 6;              //!maximum number of L_init
	iterate.n_init = 0;
	L_ini_min = 4;
	
	if(aligns.n_ali < 4)
		L_ini_min = aligns.n_ali;
	
	for(i=1; i<= n_init_max-1; i++)
	{
		iterate.n_init += 1;
		
		iterate.L_ini[iterate.n_init - 1] = aligns.n_ali / (int)pow(2.0, (iterate.n_init - 1)*1.0); //prototype double pow(double, double)
		
		if(iterate.L_ini[iterate.n_init-1] <= L_ini_min)
		{
			iterate.L_ini[iterate.n_init-1] = L_ini_min;
			return 0;
		}
	}
	iterate.n_init += 1;
	iterate.L_ini[iterate.n_init - 1] = L_ini_min;
	return 0;
}



int findMaxScore(StruInfo &stru, SeqInfo &seqs,  AlignInfo &aligns, DPara &dPara, ScoresInfo &scores, GdtInfo &gdt, IterateInfo &iterate, RmsdInfo &rmsds)
{
	int i, j, k, it, m;
	int i_init, L_init, iL_max, iL;
	int LL, ka, neq;
	double w[NMAX], r_1[3][NMAX], r_2[3][NMAX], rms;
	double u[3][3], t[3];
	int ier;

	for(i=0; i<NMAX; i++)
		w[i] = 1.0;
	
	for (i=1; i<=NMAX; i++)
	{
		for(j = 1; j<=3; j++)
		{
			r_1[j-1][i - 1] = 0.0;
			r_2[j-1][i - 1] = 0.0;
		}
	}

	//******************************************************************
	//*     find the maximum score starting from local structures superposition
	//******************************************************************
	scores.score_max = -1;              //!TM-score
	scores.score_maxsub_max = -1;       //!MaxSub-score
	scores.score10_max = -1;            //!TM-score10

	gdt.n_GDT05_max = -1;            //!number of residues<0.5
	gdt.n_GDT1_max = -1;             //!number of residues<1
	gdt.n_GDT2_max = -1;             //!number of residues<2
	gdt.n_GDT4_max = -1;             //!number of residues<4
	gdt.n_GDT8_max = -1;             //!number of residues<8
	
	for (i_init = 1; i_init <= iterate.n_init; i_init++)
	{
        L_init = iterate.L_ini[i_init-1];
		
        iL_max = aligns.n_ali - L_init + 1;
		
        for( iL=1; iL<=iL_max; iL++)       //on aligned residues, [1,nseqA]
        {
			LL = 0;
			ka = 0;
            for(i=1; i<=L_init; i++)
			{
				k = iL + i - 1;          //[1,n_ali] common aligned
				r_1[0][i-1] = stru.xa[aligns.iA[k-1]-1];
				r_1[1][i-1] = stru.ya[aligns.iA[k-1]-1];
				r_1[2][i-1] = stru.za[aligns.iA[k-1]-1];

				r_2[0][i-1] = stru.xb[aligns.iB[k-1]-1];
				r_2[1][i-1] = stru.yb[aligns.iB[k-1]-1];
				r_2[2][i-1] = stru.zb[aligns.iB[k-1]-1];

				ka = ka + 1;
				aligns.k_ali[ka-1] = k;
				LL = LL + 1;
			}
			
			u3b(w, r_1, r_2, LL, 1, rms, u, t, ier);  //u rotate r_1 to r_2
			
			if( i_init == 1 )      //global superposition
			{
				rmsds.armsd = (double)sqrt((double)(rms/LL));
				rmsds.rmsd_ali = rmsds.armsd;
			}
			
			for(j=1; j<= seqs.nseqA; j++)
			{
				stru.xt[j-1] = (double)(t[0] + u[0][0] * stru.xa[j-1] + u[0][1] * stru.ya[j-1] + u[0][2] * stru.za[j-1]);
				stru.yt[j-1] = (double)(t[1] + u[1][0] * stru.xa[j-1] + u[1][1] * stru.ya[j-1] + u[1][2] * stru.za[j-1]);
				stru.zt[j-1] = (double)(t[2] + u[2][0] * stru.xa[j-1] + u[2][1] * stru.ya[j-1] + u[2][2] * stru.za[j-1]);
			}
			dPara.d = (double)(dPara.d0_search - 1);
			
			score_fun(stru, seqs, aligns, dPara, scores, gdt);       //init, get scores, n_cut+i_ali(i) for iteration
			
			if(scores.score_max < scores.score)
			{
				scores.score_max = scores.score;
				aligns.ka0 = ka;
				for(i = 1; i<= aligns.ka0; i++)
				{
					aligns.k_ali0[i-1] = aligns.k_ali[i-1];
				}
			}
			
			if(scores.score10_max < scores.score10)
				scores.score10_max = scores.score10;
			
			if(scores.score_maxsub_max < scores.score_maxsub)
				scores.score_maxsub_max = scores.score_maxsub;
			
			if(gdt.n_GDT05_max < gdt.n_GDT05)
				gdt.n_GDT05_max = gdt.n_GDT05;
			
			if(gdt.n_GDT1_max < gdt.n_GDT1)
				gdt.n_GDT1_max = gdt.n_GDT1;
			
			if(gdt.n_GDT2_max < gdt.n_GDT2)
				gdt.n_GDT2_max = gdt.n_GDT2;
			
			if(gdt.n_GDT4_max < gdt.n_GDT4)
				gdt.n_GDT4_max = gdt.n_GDT4;
			
			if(gdt.n_GDT8_max < gdt.n_GDT8)
				gdt.n_GDT8_max = gdt.n_GDT8;
			
			
			//   iteration for extending
			dPara.d = (double)(dPara.d0_search+1);
			for(it = 1; it <= iterate.n_it; it++)
			{
				LL = 0;
				ka = 0;
				for(i=1; i<= aligns.n_cut; i++)
				{
					m = aligns.i_ali[i-1];     //[1,n_ali]
					r_1[0][i-1] = stru.xa[aligns.iA[m-1]-1];
					r_1[1][i-1] = stru.ya[aligns.iA[m-1]-1];
					r_1[2][i-1] = stru.za[aligns.iA[m-1]-1];
					r_2[0][i-1] = stru.xb[aligns.iB[m-1]-1];
					r_2[1][i-1] = stru.yb[aligns.iB[m-1]-1];
					r_2[2][i-1] = stru.zb[aligns.iB[m-1]-1];
					ka = ka+1;
					aligns.k_ali[ka-1] = m;
					LL=LL+1;
				}
				
				u3b(w, r_1, r_2, LL, 1, rms, u, t, ier);  //u rotate r_1 to r_2
				
				for(j=1; j<=seqs.nseqA; j++)
				{
					stru.xt[j-1] = (double)(t[0] + u[0][0] * stru.xa[j-1] + u[0][1] * stru.ya[j-1] + u[0][2] * stru.za[j-1]);
					stru.yt[j-1] = (double)(t[1] + u[1][0] * stru.xa[j-1] + u[1][1] * stru.ya[j-1] + u[1][2] * stru.za[j-1]);
					stru.zt[j-1] = (double)(t[2] + u[2][0] * stru.xa[j-1] + u[2][1] * stru.ya[j-1] + u[2][2] * stru.za[j-1]);
				}
				
				score_fun(stru, seqs, aligns, dPara, scores, gdt);    //get scores, n_cut+i_ali(i) for iteration
				
				if(scores.score_max < scores.score)
				{
					scores.score_max = scores.score;
					aligns.ka0 = ka;
					for(i=1; i<=ka; i++)
					{
						aligns.k_ali0[i-1] = aligns.k_ali[i-1];
					}
				}
				
				if(scores.score10_max < scores.score10)
					scores.score10_max = scores.score10;
				
				if(scores.score_maxsub_max < scores.score_maxsub)
					scores.score_maxsub_max = scores.score_maxsub;
				
				if(gdt.n_GDT05_max < gdt.n_GDT05)
					gdt.n_GDT05_max = gdt.n_GDT05;
				
				if(gdt.n_GDT1_max < gdt.n_GDT1)
					gdt.n_GDT1_max = gdt.n_GDT1;
				
				if(gdt.n_GDT2_max < gdt.n_GDT2)
					gdt.n_GDT2_max = gdt.n_GDT2;
				
				if(gdt.n_GDT4_max < gdt.n_GDT4)
					gdt.n_GDT4_max = gdt.n_GDT4;
				
				if(gdt.n_GDT8_max < gdt.n_GDT8)
					gdt.n_GDT8_max = gdt.n_GDT8;
				
				if(it == iterate.n_it)
					break;
				
				if(aligns.n_cut == ka)
				{
					neq = 0;
					for(i=1; i<=aligns.n_cut; i++)
					{
						if(aligns.i_ali[i-1] == aligns.k_ali[i-1])
							neq = neq + 1;
					}
					
					if(aligns.n_cut == neq)
						break;
				}
			}
		}
	}
	return 0;
}
		

int outputResults(StruInfo &stru, SeqInfo &seqs, AlignInfo &aligns, ScoresInfo &scores, GdtInfo &gdt, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, double plfSuperMatrix[3][4])
{
	//   output TM-scale ---------------------------->*/
	
	double score_GDT, score_GDT_HA;
	int LL, i, j;
	FILE *fp;
	//FILE *fp_test;
	
	//parameters for u3b
	double r_1[3][NMAX], r_2[3][NMAX], t[3], u[3][3], w[NMAX], rms;
	int m, ier;
	
	/*initialize weight and matrix*/
	for (i=0; i<NMAX; ++i)
		w[i] = 1.0;
	
	for (i = 1; i <= NMAX; ++i)
	{
		for(j = 1; j<=3; j++)
		{
			r_1[j-1][i - 1] = 0.0;
			r_2[j-1][i - 1] = 0.0;
		}
	}
	
	if (inputPara.m_outResult == 1)
	{	
		if ((fp = fopen(inputPara.resultFileName, "a+")) == NULL)
		{
			printf("Error: open file %s failed", inputPara.resultFileName);
			return 1;
		}
	}
	else if(inputPara.m_outResult == 0)
		fp = stderr;
	else
		fp = NULL;
		

	score_GDT = (gdt.n_GDT1_max + gdt.n_GDT2_max + gdt.n_GDT4_max + gdt.n_GDT8_max)/(double)(4*seqs.nseqB);
	score_GDT_HA = (gdt.n_GDT05_max + gdt.n_GDT1_max + gdt.n_GDT2_max + gdt.n_GDT4_max)/(double)(4*seqs.nseqB);

	if (inputPara.m_outResult > -1)
	{
		fprintf(fp,"*****************************************************************************\n");
		fprintf(fp,"*                                 TM-SCORE                                  *\n");
		fprintf(fp,"* A scoring function to assess the quality of protein structure predictions *\n");
		fprintf(fp,"* Based on statistics:                                                      *\n");
		fprintf(fp,"*       0.0 < TM-score < 0.17, Random predictions                           *\n");
		fprintf(fp,"*       0.4 < TM-score < 1.00, Meaningful predictions                       *\n");
		fprintf(fp,"* Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *\n");
		fprintf(fp,"* For comments, please email to: yzhang@ku.edu                              *\n");
		fprintf(fp,"*****************************************************************************\n");
		
		/* output TM-score*/
		fprintf(fp, "Structure1: %s Length=%4d\n", inputPara.pdb[0], seqs.nseqA);
		fprintf(fp, "Structure2: %s Length=%4d (by which all scores are normalized\n", inputPara.pdb[1], seqs.nseqB);
		
		fprintf(fp, "Number of residues in common=%4d\n", aligns.n_ali);
		fprintf(fp, "RMSD of the common residues=%8.3f\n", rmsds.rmsd_ali);
		
		fprintf(fp, "TM-score     = %6.4f (d0=%5.2f, TM10=%6.4f)\n", scores.score_max, dPara.d0, scores.score10_max);
		fprintf(fp, "MaxSub-score = %6.4f, (d0=3.50)\n", scores.score_maxsub_max);
			
		fprintf(fp, "GDT-TS-score = %6.4f %(d<1)=%6.4f, %(d<2)=%6.4f, %(d<4)=%6.4f, %(d<8)=%6.4f\n",
			score_GDT,
			gdt.n_GDT1_max/(double)(seqs.nseqB),
			gdt.n_GDT2_max/(double)(seqs.nseqB),
			gdt.n_GDT4_max/(double)(seqs.nseqB),
			gdt.n_GDT8_max/(double)(seqs.nseqB));		
		
		fprintf(fp, "GDT-HA-score = %6.4f %(d<0.5)=%6.4f, %(d<1)=%6.4f, %(d<2)=%6.4f, %(d<4)=%6.4f\n",
			score_GDT_HA,
			gdt.n_GDT05_max/(double)(seqs.nseqB),
			gdt.n_GDT1_max/(double)(seqs.nseqB),
			gdt.n_GDT2_max/(double)(seqs.nseqB),
			gdt.n_GDT4_max/(double)(seqs.nseqB));
	}
	
	//fprintf(fp_test,"%8.2f %8.2f\n", aligns.n_ali*100/(double)(seqs.nseqB), 100*(1 - aligns.n_ali/(double)(seqs.nseqB)));
	//fprintf(fp_test,"%-8.3f %-8.3f %-8.3f %-8.3f \n\n",scores.score_max, score_GDT, score_GDT_HA, rmsds.rmsd_ali);
	//fclose(fp_test);
	//   recall and output the superposition of maxiumum TM-score:
	LL=0;
	for(i=1; i<= aligns.ka0; i++)
	{
		m = aligns.k_ali0[i-1];            //!record of the best alignment
		r_1[0][i-1]=stru.xa[aligns.iA[m-1]-1];
		r_1[1][i-1]=stru.ya[aligns.iA[m-1]-1];
		r_1[2][i-1]=stru.za[aligns.iA[m-1]-1];
		r_2[0][i-1]=stru.xb[aligns.iB[m-1]-1];
		r_2[1][i-1]=stru.yb[aligns.iB[m-1]-1];
		r_2[2][i-1]=stru.zb[aligns.iB[m-1]-1];
		LL=LL+1;
	}
	
	u3b(w, r_1, r_2, LL, 1, rms, u, t, ier);  //!u rotate r_1 to r_2
	
	for(j=1; j<= seqs.nseqA; j++)
	{
		stru.xt[j-1] = (double)(t[0] + u[0][0] * stru.xa[j-1] + u[0][1] * stru.ya[j-1] + u[0][2] * stru.za[j-1]);
		stru.yt[j-1] = (double)(t[1] + u[1][0] * stru.xa[j-1] + u[1][1] * stru.ya[j-1] + u[1][2] * stru.za[j-1]);
		stru.zt[j-1] = (double)(t[2] + u[2][0] * stru.xa[j-1] + u[2][1] * stru.ya[j-1] + u[2][2] * stru.za[j-1]);
	}
	
	if (inputPara.m_outResult > -1)
	{
		fprintf(fp, "------- rotation matrix to rotate Chain-1 to Chain-2 ------\n");
		fprintf(fp, " i         t(i)          u(i,1)          u(i,2)         u(i,3)\n");
		for (i=1; i<=3; i++)
		{
			fprintf(fp, "%2d%18.10f%15.10f%15.10f%15.10f\n",
				i, t[i-1], u[i-1][0],u[i-1][1],u[i-1][2]);
		}
		
		if (inputPara.m_outResult == 1)
			fclose(fp);
	}
	
	for (i=0; i<3; i++)
	{
		plfSuperMatrix[i][0] =  u[i][0];
		plfSuperMatrix[i][1] =  u[i][1];
		plfSuperMatrix[i][2] =  u[i][2];
		plfSuperMatrix[i][3] =  t[i];
	}
	
	
	// rmsd in superposed regions
	dPara.d = (double)(dPara.d_output);                //for output
	score_fun(stru, seqs, aligns, dPara, scores, gdt);                              //give i_ali(i), score_max=score now
	
	LL=0;
	for(i=1; i<= aligns.n_cut; i++)
	{
		m = aligns.i_ali[i-1];             //[1,nseqA]
		r_1[0][i-1] = stru.xa[aligns.iA[m-1]-1];
		r_1[1][i-1] = stru.ya[aligns.iA[m-1]-1];
		r_1[2][i-1] = stru.za[aligns.iA[m-1]-1];
		r_2[0][i-1] = stru.xb[aligns.iB[m-1]-1];
		r_2[1][i-1] = stru.yb[aligns.iB[m-1]-1];
		r_2[2][i-1] = stru.zb[aligns.iB[m-1]-1];
		LL=LL+1;
	}
	
	u3b(w, r_1, r_2, LL, 0, rms, u, t, ier);
	rmsds.armsd = (double)sqrt((double)(rms/LL));
	rmsds.rmsd = rmsds.armsd;	
	return 0;
}


/////////////////////////////////////////////////////////////////////////////
// Member Function: outputRotatedChain
// Purpose:  	
// Parameters:
// Returns:
// Comments:
/////////////////////////////////////////////////////////////////////////////
int outputRotatedChain(StruInfo &stru,  AlignInfo &aligns, ScoresInfo &scores, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, SeqInfo &seqs)
{
//   output rotated chain1 + chain2----->
      	
	int i;
	FILE* fp;
	char tempData[100];
	
	if (inputPara.m_out != 1)
		return 1;
	
	if ((fp = fopen(inputPara.outname, "w")) == NULL)  //pdb1.aln + pdb2.aln
	{
		printf("\nError: open file %s failed", inputPara.outname);
		return 1;
	}

 	fprintf(fp, "load inline\n");
	fprintf(fp, "select atomno<1000\n");
	fprintf(fp, "wireframe .45\n");
	fprintf(fp, "select none\n");
	fprintf(fp, "select atomno>1000\n");
	fprintf(fp, "wireframe .15\n");
	fprintf(fp, "color white\n");

    for (i = 1; i <= aligns.n_cut; ++i)
	{
		fprintf(fp, "select %4d\n", seqs.nresA[aligns.iA[aligns.i_ali[i - 1] - 1] - 1]);
		fprintf(fp, "color red\n");
	}

	fprintf(fp, "select all\n");
	fprintf(fp, "exit\n");
	fprintf(fp, "REMARK  RMSD of the common residues=%8.3f\n",  rmsds.rmsd_ali);
	fprintf(fp, "REMARK TM-score=%6.4f (d0= %5.2f)\n", scores.score_max, dPara.d0);
	
	for (i = 1; i <= seqs.nseqA; ++i)
	{
     	fprintf(fp, "ATOM  %5d  CA  ", seqs.nresA[i - 1]);
		sprintf(tempData, "%.*s", 3, seqs.seqA + (i - 1) * 3);
		fprintf(fp, "%s%6d    %8.3f%8.3f%8.3f\n",
			tempData, seqs.nresA[i - 1], stru.xt[i - 1], stru.yt[i - 1], stru.zt[i - 1]);
	}
	fprintf(fp, "TER\n");

	for (i = 2; i <= seqs.nseqA; ++i)
	{
		fprintf(fp, "CONECT%5d%5d\n", seqs.nresA[i - 2], seqs.nresA[i - 1]);
	}

	for (i = 1; i <= seqs.nseqB; ++i)
	{
		fprintf(fp, "ATOM  %5d  CA  ", seqs.nresB[i - 1] + 2000);
		sprintf(tempData, "%.*s", 3, seqs.seqB + (i - 1) * 3);
		fprintf(fp, "%s%6d    %8.3f%8.3f%8.3f\n",
			tempData, seqs.nresB[i - 1], stru.xb[i - 1], stru.yb[i - 1], stru.zb[i - 1]);
	}
	fprintf(fp, "TER\n");

	for (i = 2; i <= seqs.nseqB; ++i)
	{
		fprintf(fp, "CONECT%5d%5d\n", seqs.nresB[i - 2] + 2000, seqs.nresB[i - 1] + 2000);
	}


	for (i = 2; i <= seqs.nseqA; ++i)
		fprintf(fp, "CONECT%5d%5d\n", seqs.nresA[i-2], seqs.nresA[i-1]);

	
    for (i = 2; i <= seqs.nseqB; ++i)
	{
		fprintf(fp, "ATOM  %5d  CA  ", seqs.nresB[i - 1]+2000);
		sprintf(tempData, "%.*s", 3, seqs.seqB + (i - 1) * 3);
		fprintf(fp, "%s%6d    %8.3f%8.3f%8.3f\n",
			tempData, seqs.nresB[i - 1], stru.xb[i - 1], stru.yb[i - 1], stru.zb[i - 1]);
	}

	fprintf(fp, "TER\n");
	for (i = 2; i <= seqs.nseqB; ++i)
		fprintf(fp, "CONECT%5d%5d\n", seqs.nresB[i - 2]+2000, 2000+seqs.nresB[i-1]);

	fclose(fp);

	return 0;
}


/////////////////////////////////////////////////////////////////////////////
// Member Function: outputAlignedSequences
// Purpose:  	
// Parameters:
// Returns:
// Comments:
/////////////////////////////////////////////////////////////////////////////
int outputAlignedSequences(StruInfo &stru,  AlignInfo &aligns, InputPara &inputPara, DPara &dPara, RmsdInfo &rmsds, SeqInfo &seqs)
{
	//record aligned residues by i=[1,nseqA], for sequenceM()------------>	
	int i, j, k;
	char sequenceA[1*NMAX], sequenceB[1*NMAX], sequenceM[1*NMAX];
	int iq[NMAX];
	double r1, r2, r3, dis;
	FILE *fp;
	
	if (inputPara.m_outResult !=1 )
		return 1;
	
	if ((fp = fopen(inputPara.resultFileName, "a+")) == NULL)
	{
		printf("Error: append file %s failed\n", inputPara.resultFileName);
		return 1;
	}
	
	
	for (i = 1; i <= seqs.nseqA; ++i)
		iq[i - 1] = 0;
	
	for (i = 1; i <= aligns.n_cut; ++i)
	{
		j = aligns.iA[aligns.i_ali[i-1]-1];         //[1,nseqA]
		k = aligns.iB[aligns.i_ali[i-1]-1];  		 //[1,nseqB]
		r1 = stru.xt[j - 1] - stru.xb[k - 1];
		r2 = stru.yt[j - 1] - stru.yb[k - 1];
		r3 = stru.zt[j - 1] - stru.zb[k - 1];
		
		dis = (double)sqrt(r1*r1 + r2*r2 + r3*r3);
		
		if(dis < dPara.d_output)
            iq[j - 1] = 1;
		
	}
	
	//output aligned sequences
	k = 0;
	i = 1;
	j = 1;
	for(;;)
	{
		if(i > seqs.nseqA && j > seqs.nseqB)
			break;
		if(i > seqs.nseqA && j <= seqs.nseqB)
		{
			k = k+1;
			sequenceA[k-1] = '-';
			sequenceB[k-1] = seqs.seq1B[j-1];
			sequenceM[k-1] = ' ';
			j = j+1;
			continue;
		}
		
		if(i <= seqs.nseqA && j > seqs.nseqB)
		{
			k = k+1;
			sequenceA[k-1] = seqs.seq1A[i-1];
			sequenceB[k-1] = '-';
			sequenceM[k-1] = ' ';
			i = i+1;
			continue;
		}
		
		if(seqs.nresA[i-1] == seqs.nresB[j-1])
		{
			k = k+1;
			sequenceA[k-1] = seqs.seq1A[i-1];
			sequenceB[k-1] = seqs.seq1B[j-1];
			if(iq[i-1] == 1)
				sequenceM[k-1] = ':';
			else
				sequenceM[k-1]= ' ';
			
			i = i+1;
			j = j+1;
			continue;
		}
		else if( seqs.nresA[i-1] < seqs.nresB[j-1])
		{
			k = k+1;
			sequenceA[k-1] = seqs.seq1A[i-1];
			sequenceB[k-1] = '-';
			sequenceM[k-1] = ' ';
			i = i+1;
			continue;
		}
		else if( seqs.nresB[j-1] < seqs.nresA[i-1])
		{
			k = k+1;
			sequenceA[k-1] = '-';
			sequenceB[k-1] = seqs.seq1B[j-1];
			sequenceM[k-1] = ' ';
			j = j+1;
			continue;
		}
	}
	
   	fprintf(fp, "Superstation in the TM-score: Length (d<%3.1f )=%d RMSD=%6.2f\n",
		dPara.d_output, aligns.n_cut,  rmsds.rmsd);
	
	fprintf(fp, "(\":\" denotes the residue pairs of distance<%4.1f Angstrom\n", dPara.d_output);
	
    sequenceA[k] = '\0';
	sequenceM[k] = '\0';
	sequenceB[k] = '\0';

	fprintf(fp, "%s\n", sequenceA);	
	fprintf(fp, "%s\n", sequenceM);	
	fprintf(fp, "%s\n", sequenceB);
	
	for (i = 1; i <= k; ++i)
		fprintf(fp, "%d", i%10);
	fprintf(fp, "\n\n");
	
	fclose(fp);
	return 0;
}


/////////////////////////////////////////////////////////////////////////////
// Member Function: score_fun
// Purpose: 1, collect those residues with dis<d  	
//			2, calculate score_GDT, score_maxsub, score_TM
// Parameters:
// Returns:
// Comments:
/////////////////////////////////////////////////////////////////////////////
int score_fun(StruInfo &stru, SeqInfo &seqs, AlignInfo &aligns, DPara &dPara, ScoresInfo &scores, GdtInfo &gdt)
{
	double d_tmp;
	double score_maxsub_sum, score_sum, score_sum10;
	int i, j, k;
	double r1, r2, r3, dis;

	d_tmp = dPara.d;	
	do
	{
		aligns.n_cut = 0;                //number of residue-pairs dis<d, for iteration
		gdt.n_GDT05 = 0;                 //for GDT-score, # of dis<0.5
		gdt.n_GDT1 = 0;                  //for GDT-score, # of dis<1
		gdt.n_GDT2 = 0;                  //for GDT-score, # of dis<2
		gdt.n_GDT4 = 0;                  //for GDT-score, # of dis<4
		gdt.n_GDT8 = 0;                  //for GDT-score, # of dis<8
		
		score_maxsub_sum = 0;        //Maxsub-score
		score_sum = 0;               //TMscore
		score_sum10 = 0;             //TMscore10
		
		for(k = 1; k <= aligns.n_ali; ++k)
		{
			i = aligns.iA[k-1];                //[1,nseqA] reoder number of structureA
			j = aligns.iB[k-1];                //[1,nseqB]
			r1 = stru.xt[i - 1] - stru.xb[j - 1];
			r2 = stru.yt[i - 1] - stru.yb[j - 1];
			r3 = stru.zt[i - 1] - stru.zb[j - 1];
			dis = sqrt(r1 * r1 + r2 * r2 + r3 * r3);
			//for iteration:
			if(dis < d_tmp)
			{
				aligns.n_cut = aligns.n_cut + 1;
				aligns.i_ali[aligns.n_cut - 1] = k;      //[1,n_ali], mark the residue-pairs in dis<d
			}
			//   for GDT-score:
			if(dis <= 8)
			{
				gdt.n_GDT8 = gdt.n_GDT8 + 1;
				if(dis <= 4)
				{
					gdt.n_GDT4 = gdt.n_GDT4 + 1;
					if(dis <= 2)
					{
						gdt.n_GDT2 = gdt.n_GDT2 + 1;
						if(dis <= 1)
						{
							gdt.n_GDT1 = gdt.n_GDT1 + 1;
							if(dis <= 0.5)
								gdt.n_GDT05 = gdt.n_GDT05 + 1;
						}
					}
				}
			}
			//   for MAXsub-score:
			if(dis < 3.5)
			{
				r1 = dis/3.5;
				score_maxsub_sum += 1/(r1 * r1 + 1);			
			}
			
			//for TM-score:
			r1 = dis/dPara.d0;
			score_sum +=  1/(1+r1*r1);
			
			//for TM-score10:
			if(dis < 10)
			{
				r1 = dis/dPara.d0;
				score_sum10 += 1/(1+r1*r1);
			}
		}
		d_tmp += 0.5;
		
	}while(aligns.n_cut < 3 && aligns.n_ali > 3);
	
	scores.score_maxsub = score_maxsub_sum /(double)(seqs.nseqB); //MAXsub-score
	scores.score = score_sum/(double)(seqs.nseqB);                //TM-score
	scores.score10 = score_sum10/(double)(seqs.nseqB);            //TM-score10
	
	return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Member Function: u3b
// Purpose:
// Parameters:
// Returns:
// Comments:
/////////////////////////////////////////////////////////////////////////////
/*************************************************************************
*                                 u3b                                   *
*  Calculate sum of (r_d-r_m)^2                                         *
*  w    - w(m) is weight for atom pair  c m                 (given)     *
*  x    - x(i,m) are coordinates of atom c m in set x       (given)     *
*  y    - y(i,m) are coordinates of atom c m in set y       (given)     *
*  n    - n is number of atom pairs                         (given)     *
*  mode  - 0:calculate rms only                             (given)     *
*          1:calculate rms,u,t                           (takes longer) *
*  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)    *
*  u    - u(i,j) is   rotation  matrix for best superposition  (result) *
*  t    - t(i)   is translation vector for best superposition  (result) *
*  ier  - 0: a unique optimal superposition has been determined(result) *
*       -1: superposition is not unique but optimal                     *
*       -2: no result obtained because of negative weights w            *
*           or all weights equal to zero.                               *
*************************************************************************/
//subroutine u3b(w, x, y, n, mode, rms, u, t, ier)

int u3b(double w[NMAX], double x[3][NMAX], double y[3][NMAX], int n, int mode, double &rms, double u[3][3], double t[3], int &ier)
{
	int i, j, k, l, m1, m;
		
	double r[3][3], xc[3], yc[3], wc, a[3][3], b[3][3], e0;
	double d, h, g, p;		
	double e[3], rr[6], ss[6];		
	double sqrt3 = 1.73205080756888;
	double tol = 0.01;
	double zero = 0.0;
	double one = 1.0, two = 2.0, three = 3.0;
	int ip[9] = { 1, 2, 4, 2, 3, 5, 4, 5, 6};
	int ip2312[4] = { 2, 3, 1, 2 };
		
	double cof, det, cth, sth, spur, sigma, sqrth;
	double d1, d2;
		
	// "rms.for"   //zjf_check
		
	wc = zero;		
	rms = 0.0;
	e0 = zero;
		
	for (i = 1; i <= 3; i++) 		
	{
		xc[i-1] = zero;
		yc[i-1] = zero;
		t[i-1] = 0.0;			
		
		for (j = 1; j <= 3; j++)
		{		
			d = zero;			
			if (i == j)			
				d = one;				
			u[i-1][j-1] = d;			
			a[i-1][j-1] = d;			
			r[i-1][j-1] = zero;			
		}		
	}	
	
	ier = -1;	
	/**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y*/		  	
	// 170 "rms.for"	
	if (n < 1)	
		return 0;
	
	// 172 "rms.for"
	ier = -2;	
	for (m = 1; m <= n; m++) 	
	{
		if (w[m-1] < 0.0) 		
			return 0;			
		
		wc += w[m-1];		
		
		for (i = 1; i <= 3; i++)              //zjf_check		
		{
			xc[i-1] += w[m-1] * x[i-1][m-1];
			yc[i-1] += w[m-1] * y[i-1][m-1];
		}
	}
	
	if (wc <= zero)
		return 0;
	
	for (i = 1; i <= 3; ++i) 		
	{
		xc[i - 1] /= wc;
		/* **** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X */
		/* 182 "rms.for" */
		/* L3: */
		yc[i - 1] /= wc;
		
	}
		
	
	// 184 "rms.for"	
	for (m = 1; m <= n; ++m)
	{
		for (i = 1; i <= 3; ++i)
		{
			/* Computing 2nd power */
			d1 = x[i-1][m-1] - xc[i - 1];
			/* Computing 2nd power */
			d2 = y[i-1][m-1] - yc[i - 1];
			e0 += w[m] * (d1 * d1 + d2 * d2);
			/* 187 "rms.for" */
			d = w[m-1] * y[i-1][m-1] - yc[i - 1];
			for (j = 1; j <= 3; j++)
			{
				/* **** CALCULATE DETERMINANT OF R(I,J) */
				/* 189 "rms.for" */
				/* L4: */
				r[i-1][j-1] += d * (x[j-1][m-1] - xc[j - 1]);
			}
		}
	}
	
	/* 191 "rms.for" */
	det = r[0][0] * ((r[1][1] * r[2][2]) - (r[1][2] * r[2][1])) - 	
		  r[0][1] * ((r[1][0] * r[2][2]) - (r[1][2] * r[2][0])) +
		  r[0][2] * ((r[1][0] * r[2][1]) - (r[1][1] * r[2][0]));
	
	/* **** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R */
	/* 194 "rms.for" */
	sigma = det;
	/* 196 "rms.for" */
	m = 0;
	
	for (j = 1; j <= 3; ++j)
	{
		for (i = 1; i <= j; ++i)
		{
			m++;
			/* ***************** EIGENVALUES ***************************************** */
			/* **** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0 */
			/* 200 "rms.for" */
			/* L5: */
			rr[m-1] = (r[0][i-1] * r[0][j-1]) + (r[1][i-1] * r[1][j-1]) + (r[2][i-1] * r[2][j-1]);
		}
	}
	// 203 "rms.for"
	spur = ((rr[0] + rr[2]) + rr[5])/three;
	cof = (rr[2] * rr[5] - rr[4] * rr[4] + rr[0] * rr[5] - rr[3] * rr[3] + rr[0] * rr[2] - rr[1] * rr[1])/three;
	// 205 "rms.for"
	det = det * det;
	for (i = 1; i <= 3; ++i)
	{
		/* L6: */
		e[i - 1] = spur;
	}
	
	/* **** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR */
	/* 208 "rms.for" */
	if (spur < zero)
		goto L40;
	
	// 210 "rms.for"
	d = spur * spur;
	h = d - cof;
	
	/* **** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER */
	/* 212 "rms.for" */
	g = (spur * cof - det)/two - spur * h;
	/* 214 "rms.for" */
	if (h <= zero)
	{
		goto L8;
	}
	
	sqrth = sqrt(h);
	d = h * h * h - g * g;
	
	if (d < zero)
		d = zero;
	
	d = atan2(sqrt(d), -g)/three;
	cth = sqrth * cos(d);
	sth = sqrth * sqrt3 * sin(d);
	
	e[0] = spur + cth + cth;
	e[1] = spur - cth + sth;
	e[2] = spur - cth - sth;
	
	/* .....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS */
	/* 224 "rms.for" */
	if (mode !=0)
		goto L10;
	else
		goto L50;
	
	/* **************** EIGENVECTORS ***************************************** */
	/* 226 "rms.for" */
L8:
	if (mode !=0)
		goto L30;
	else
		goto L50;
	/* 228 "rms.for" */
L10:	
	for(l = 1; l<=3;  l+=2)
	{
		d = e[l-1];		
		ss[0] = (d - rr[2]) * (d - rr[5]) - (rr[4] * rr[4]);
		ss[1] = (d - rr[5]) * rr[1] + (rr[3] * rr[4]);
		ss[2] = (d - rr[0]) * (d - rr[5]) - (rr[3] * rr[3]);
		ss[3] = (d - rr[2]) * rr[3] + (rr[1] * rr[4]);
		ss[4] = (d - rr[0]) * rr[4] + (rr[1] * rr[3]);
		ss[5] = (d - rr[0]) * (d - rr[2]) - (rr[1] * rr[1]);
		j = 1;
		if (fabs(ss[0]) >= fabs(ss[2]))
			goto L12;
		
		j = 2;
		if (fabs(ss[2]) >= fabs(ss[5]))
			goto L13;
		
L11:
		j = 3;
		goto L13;
		
L12:
		if (fabs(ss[0]) < fabs(ss[5]))
			goto L11;
		
L13:
		d = zero;
		j = 3 * (j - 1);
		for(i = 1;i<=3; i++)
		{
			k = ip[i+j-1];
			a[i-1][l-1] = ss[k-1];
			d += ss[k-1] * ss[k-1];
		}
		
		if (d > zero)
			d = one/sqrt(d);
		
		for(i = 1;i<=3; i++)
			a[i-1][l-1] = a[i-1][l-1] * d;
	}
		
	d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
	m1 = 3;
	m = 1;
	
	if ((e[0] - e[1]) > (e[1] - e[2]))
		goto L16;
	
	m1 = 1;
	m = 3;

L16:
	p = zero;
	for(i = 1; i<= 3; i++)
	{
		a[i-1][m1-1] -= d * a[i-1][m-1];
		p += a[i-1][m1-1]*a[i-1][m1-1];
	}
	
	if (p < tol)
		goto L19;
	
	p = one/sqrt(p);
	
	for(i = 1; i<= 3; i++)
		a[i-1][m1-1] *= p;
	
	goto L21;

L19:
	p = one;
	for(i = 1; i<= 3; i++)
	{
		if (p < fabs(a[i-1][m-1]) )
			continue;

		p = fabs(a[i-1][m-1]);
		j = i;		
	}
	
	k = ip2312[j-1];
	l = ip2312[j];
	p = sqrt(a[k-1][m-1]*a[k-1][m-1] + a[l-1][m-1]*a[l-1][m-1]);
	if (p <= tol)
		goto L40;
	
	a[j-1][m1-1] = zero;
	a[k-1][m1-1] = - (a[l-1][m-1]/p);
	a[l-1][m1-1] = a[k-1][m-1]/p;
	
L21:
	a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
	a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];	
	/* ****************** ROTATION MATRIX ************************************ */
	/* 282 "rms.for" */	
	a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
	
	// 284 "rms.for"
L30:
	for (l = 1; l <= 2; l++)
	{
		d = zero;
		
		for (i = 1; i <= 3; ++i)
		{
			b[i-1][l-1] = r[i-1][0] * a[0][l-1] + r[i-1][1] * a[1][l-1] + r[i-1][2] * a[2][l-1];
			// 288 "rms.for"
			d += b[i-1][l-1]*b[i-1][l-1];
		}
		
		if (d > zero)
			d = one/sqrt(d);
		
		for (i = 1; i <= 3; i++)
			b[i-1][l-1] = b[i-1][l-1] * d;
	}

	d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
	p = zero;
	for (i = 1; i <= 3; i++)
	{
		b[i-1][1] -= d * b[i-1][0];
		p += b[i-1][1]*b[i-1][1];
	}

	if (p < tol)
		goto L35;
	
	p = one/sqrt(p);
	for (i = 1; i <= 3; i++)
		b[i-1][1] *= p;
	
	goto L37;
	
L35:
	p = one;
	for (i = 1; i <= 3; i++)
	{
		if (p < fabs(b[i-1][0]))
			continue;

		p = fabs(b[i-1][0]);
		j = i;
	}
	
	k = ip2312[j-1];
	l = ip2312[j];
	p = sqrt( b[k-1][0]*b[k-1][0] + b[l-1][0] *b[l-1][0]);
	
	if (p < tol)
		goto L40;
	
	b[j-1][1] = zero;
	b[k-1][1] = - (b[l-1][0] / p);
	b[l-1][1] = b[k-1][0] / p;
	
L37:
	b[0][2] = (b[1][0] * b[2][1]) - (b[1][1] * b[2][0]);
	b[1][2] = (b[2][0] * b[0][1]) - (b[2][1] * b[0][0]);
	b[2][2] = (b[0][0] * b[1][1]) - (b[0][1] * b[1][0]);
	
	for (i = 1; i <= 3; i++)
	{
		for (j = 1; j <= 3; j++)
		{
			/* ****************** TRANSLATION VECTOR ********************************* */
			/* 320 "rms.for" */
			/* L39: */
			u[i-1][j-1] = ((b[i-1][0] * a[j-1][0]) + (b[i-1][1] * a[j-1][1])) + (b[i-1][2] * a[j-1][2]);
		}
	}
L40:
	for (i = 1; i <= 3; i++)
	{
		/* ********************** RMS ERROR ************************************** */
		/* 323 "rms.for" */
		/* L41: */
		t[i-1] = yc[i-1] - u[i-1][0] * xc[0] - u[i-1][1] * xc[1] - u[i-1][2] * xc[2];
	}
	
L50:
	for (i = 1; i <= 3; i++)
	{
		if (e[i-1] < zero)
			e[i-1] = zero;
		e[i-1] = sqrt(e[i-1]);
	}
	
	ier = 0;
	if (e[1] < e[0] * 1e-5)
		ier = -1;
	
	d = e[2];
	if (sigma >= 0.0)
		goto L52;
	
	d = - d;
	if ((e[1] - e[2]) <= (e[0] * 1e-5))
		ier = -1;
	
L52:
	d = (d + e[1]) + e[0];
	rms = (e0 - d) - d;
	if (rms < 0.0)
		rms = 0.0;
	return 0;
}



/////////////////////////////////////////////////////////////////////////////
// Member Function: readFirstCaFile
// Purpose:  	
// Parameters:
// Returns:
// Comments:
/////////////////////////////////////////////////////////////////////////////
int readFirstCaFile(char* fileName, StruInfo &stru,  SeqInfo &seqs)
{
	FILE *fp;
	int i, j;
	char s[100], du[100], tempData[100];

	
	if ((fp = fopen(fileName, "r")) == NULL)
	{		
		printf("\nError: File %s does not exist when read first CA file\n", fileName);		
		return -1;		
	}	
	
	i = 0;	
	fgets(s, 100, fp);
	
	while (!feof(fp))		
	{		
		if (strncmp(s, "TER", 3) == 0)			
		{			
			break;			
		}
		
		if (strncmp(s, "ATOM", 4) == 0)			
		{			
			if ((strncmp(s+12, "CA  ", 4) == 0) ||				
				(strncmp(s+12, " CA ", 4) == 0) ||				
				(strncmp(s+12, "  CA", 4) == 0))				
			{				
				if (s[16] == ' ' || s[16] == 'A')					
				{					
					++i;					
					// Read from s with FMT_103:format(A17,A3,A2,i4,A4,3F8.3)
					
					sprintf(du, "%.*s", 17, s);					
					sprintf(seqs.seqA + (i -1)*3, "%.*s", 3, s+17);  					
					sprintf(du, "%.*s", 2, s+20);				
					
					sprintf(tempData, "%.*s", 4, s+22);					
					seqs.nresA[i-1] = atoi(tempData);					
					
					sprintf(du, "%.*s", 4, s+26, 4);					
					sprintf(tempData, "%.*s", 8, s+30);					
					stru.xa[i-1] = (double)atof(tempData);
					
					sprintf(tempData, "%.*s", 8, s+30+8);					
					stru.ya[i-1] = (double)atof(tempData);	
					
					sprintf(tempData, "%.*s", 8, s+30+8+8);						
					stru.za[i-1] = (double)atof(tempData);					
					
					for (j = -1; j <= 20; ++j) 						
					{						
						if (strncmp(seqs.seqA + (i - 1)*3, aa + (j + 1)*3, 3) == 0)							
						{						
							*(unsigned char *)&seqs.seq1A[i - 1] = *(unsigned char *)&slc[j + 1];							
						}						
					}						
				}				
			}			
		}		
		fgets(s, 100, fp);		
	} //while
	
	seqs.nseqA = i;	
	fclose(fp);	
	return 1;	
}



/////////////////////////////////////////////////////////////////////////////
// Member Function: readSecondCaFile
// Purpose:  	
// Parameters:
// Returns:
// Comments:
/////////////////////////////////////////////////////////////////////////////
int readSecondCaFile(char* fileName, StruInfo &stru,  SeqInfo &seqs)
{	

	FILE *fp;	
	int i, j;	
	char s[100], du[100], tempData[100];	
	
	if ((fp = fopen(fileName, "r")) == NULL)		
	{		
		printf("\nError: File %s does not exist when read second CA file\n", fileName);		
		return -1;		
	}	
	
	i = 0;	
	fgets(s, 100, fp);
	
	while (!feof(fp))		
	{		
		if (strncmp(s, "TER", 3) == 0)			
		{			
			break;			
		}		
		
		if (strncmp(s, "ATOM", 4) == 0)			
		{			
			if ((strncmp(s+12, "CA  ", 4) == 0) ||				
				(strncmp(s+12, " CA ", 4) == 0) ||				
				(strncmp(s+12, "  CA", 4) == 0))				
			{				
				if (s[16] == ' ' || s[16] == 'A')					
				{					
					++i;					
					sprintf(du, "%.*s", 17, s);					
					sprintf(seqs.seqB + (i -1)*3, "%.*s", 3, s+17);  					
					sprintf(du, "%.*s", 2, s+17+3);					
					
					sprintf(tempData, "%.*s", 4, s+17+3+2);					
					seqs.nresB[i-1] = atoi(tempData);					
					sprintf(du, "%.*s", 4, s+26, 4);					
					sprintf(tempData, "%.*s", 8, s+30);						
					stru.xb[i-1] = (double)atof(tempData);
					
					sprintf(tempData, "%.*s", 8, s+30+8);					
					stru.yb[i-1] = (double)atof(tempData);
					
					sprintf(tempData, "%.*s", 8, s+30+8+8);					
					stru.zb[i-1] = (double)atof(tempData);					
					
					for (j = -1; j <= 20; ++j) 						
					{						
						if (strncmp(seqs.seqB + (i - 1)*3, aa + (j + 1)*3, 3) == 0)							
						{							
							*(unsigned char *)&seqs.seq1B[i - 1] = *(unsigned char *)&slc[j + 1];							
						}						
					}						
				}				
			}			
		}		
		fgets(s, 100, fp);		
	} //while	
	seqs.nseqB = i;	
	fclose(fp);
	return 1;	
}



/////////////////////////////////////////////////////////////////////////////
// Member Function: pickupAlignedRes
// Purpose:  	
// Parameters:
// Returns:
// Comments:
/////////////////////////////////////////////////////////////////////////////
int 
pickupAlignedRes( AlignInfo &aligns, SeqInfo &seqs, StruInfo &stru)
{
	int i, j, k;	
	k = 0;
	for (i = 1; i <= seqs.nseqA; ++i)
	{
		for (j = 1; j <= seqs.nseqB; ++j)
		{
			if (seqs.nresA[i-1] == seqs.nresB[j-1])
			{
				++k;
				aligns.iA[k-1] = i;
				aligns.iB[k-1] = j;
				break;
			}
		}
	}	
	aligns.n_ali = k; /* number of aligned residues */
	
					 /*	k = (seqs.nseqA > seqs.nseqB)? seqs.nseqB : seqs.nseqA;
					 aligns.n_ali = k;
					 for(i = 1; i<=k; i++)
					 {
					 aligns.ia[i-1] = i;
					 aligns.ib[i-1] = i;
					 }
	*/
	if (aligns.n_ali < 1)
	{
		printf("There is no common residues in the input structures.\n");
		return -1;
	}
	return 1;
	
	
}





/*********************************************************************************
*          NAME:  P2D_genPred2D()
*   DESCRIPTION:  The scheme to evoke the prediction routines.
*         INPUT:  None
*        OUTPUT:  None
*        RETURN:  None
**********************************************************************************/
void 
transSeq2Seq3(char *ss1, char *ss3, int nSeqLen)
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

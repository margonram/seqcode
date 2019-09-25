/*************************************************************************
*                                                                        *
*   Module: seqcode.h                                                    *
*                                                                        *
*   Definitions, data structures types and imported headers              *
*                                                                        *
*   This file is part of the seqcode 1.0 distribution                    *
*                                                                        *
*     Copyright (C) 2020 - Enrique BLANCO GARCIA                         *
*                                                                        *
*  This program is free software; you can redistribute it and/or modify  *
*  it under the terms of the GNU General Public License as published by  *
*  the Free Software Foundation; either version 2 of the License, or     *
*  (at your option) any later version.                                   *
*                                                                        *
*  This program is distributed in the hope that it will be useful,       *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*  GNU General Public License for more details.                          *
*                                                                        *
*  You should have received a copy of the GNU General Public License     *
*  along with this program; if not, write to the Free Software           *
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             *
*************************************************************************/     

/* Required libraries */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

/*************************************************************************
A. DEFINITIONS
*************************************************************************/

/* Maximum chars per input line             */
#define MAXLINE 10000             

/* Maximum length for strings (mess)        */
#define MAXSTRING 1000            

/* Maximum length for gene symbol names     */
#define MAXGENENAME 100            

/* Maximum length for chromosome names      */
#define MAXCHRNAME 50

/* Maximum length for R color names         */
#define MAXCOLORNAME 50

/* Maximum length for R RGB codes (string)  */
#define MAXRGBCODE 15

/* Maximum number of R colors               */
#define MAXCOLORS 700

/* Maximum value to colorize the maps       */
#define MAXCOLORIZE 100

/* Maximum number of distinct chromosomes   */
#define MAXCHRS 100

/* Maximum number of distinct genes         */
#define MAXGENES 50000

/* Maximum length of filenames              */
#define FILENAMELENGTH 200       

/* Maximum number of usage messages         */
#define MAXUSAGE 10

/* The name of the games                    */
#define VERSION             "SeqCode_v1.0"  
#define BUILDCHIPPROFILE    "buildChIPprofile"  
#define COMBINECHIPPROFILES "combineChIPprofiles"
#define COMBINETSSMAPS      "combineTSSmaps"
#define COMBINETSSPLOTS     "combineTSSplots"  
#define COMPUTEMAXSIGNAL    "computeMaxSignal"  
#define FINDPEAKS           "findPeaks"  
#define GENOMEDISTRIBUTION  "genomeDistribution"  
#define MATCHPEAKS          "matchpeaks"  
#define MATCHPEAKSGENES     "matchpeaksgenes"
#define PROCESSMACS         "processmacs"
#define PRODUCEGENEMAPS     "produceGENEmaps"  
#define PRODUCEGENEPLOTS    "produceGENEplots"  
#define PRODUCEPEAKMAPS     "producePEAKmaps"
#define PRODUCEPEAKPLOTS    "producePEAKplots"
#define PRODUCETESMAPS      "produceTESmaps"
#define PRODUCETESPLOTS     "produceTESplots"
#define PRODUCETSSMAPS      "produceTSSmaps"  
#define PRODUCETSSPLOTS     "produceTSSplots"  
#define RECOVERCHIPLEVELS   "recoverChIPLevels"
#define SCOREPHASTCONS      "scorePhastCons"

/* URLs of SeqCode                          */
#define SOURCECODE          "https://github.com/eblancoga/seqcode"
#define HOMEPAGE            "http://ldicrocelab.crg.es"

/* Authorship                               */
#define AUTHORSHIP          "Enrique Blanco"

/* Dictionary definitions (hash)            */
#define MAXENTRY 97              
#define MAXTYPE 50               
#define MAXINFO 100
#define NOTFOUND -1

/* Tab space within fields                  */
#define TAB "\t"

/* New line character                       */
#define NL "\n"

/* White space character                    */
#define SPACE " "

/* Equal character                          */
#define sEQ "="

/* Comma between exon coordinates (refgene) */
#define CSV ","

/* Gene with no predefined transcript       */
#define NMNULL "NMNULL"

/* File size constant                       */
#define MEGABYTE 1048576

/* TMP file (BG)                            */
#define TMPFILE "/tmp/fileSeqCode"

/* UCSC headers                             */
#define TRACK "track"

/* Time constants                           */
#define MINUTE 60
#define YEAR 2020

/* Message (info)/X chars for big files     */
#define MESSAGE_FREQ 100000

/* Message (info)/X chars for small files   */
#define MESSAGE_FREQ_SMALL 1000      

/* Max number of transcripts/chr (refGene)  */
#define MAXTRANSCRIPTS 10000

/* Max number of peaks/chr                  */
#define MAXPEAKS 500000

/* DNA strands                              */
#define cFORWARD '+'
#define cREVERSE '-'

/* Boolean values                           */
#define TRUE 1
#define FALSE 0

/* Genome regions (marks: basic)            */
#define cPROMOTER 'P'
#define cINTRAGENIC 'G'

/* Genome regions (strings: basic)          */
#define sPROMOTER "PROMOTER "
#define sINTRAGENIC "INTRAGENIC "
#define sINTERGENIC "INTERGENIC "

/* Genome regions (marks: detailed)         */
#define cPROXIMAL 'X'
#define cDISTAL 'D'
#define c5UTR '5'
#define c3UTR '3'
#define cCDS 'S'
#define cINTRON 'I'

/* Genome regions (strings: detailed)       */
#define sPROXIMAL "PROXIMAL "
#define sDISTAL "DISTAL "
#define s5UTR "5UTR "
#define s3UTR "3UTR "
#define sCDS "CDS "
#define sINTRON "INTRON "

/* Exon mark (PhastCons exon filtering)     */
#define cEXON 'E'

/* Promoter limits for GenomeDistribution   */
#define PROMOTERSIZE 2500
#define PROXIMALPOINT 500

/* Non-defined scores (BED 4th column)      */
#define NOSCORE  0.0
#define sNOSCORE "0.0"

/* Range of values for plotTSS profiles     */
#define MINLENGTH 500
#define MAXLENGTH 5000

/* Window (bin) resolution for each program */                 
#define WINDOWRES_GENOME 10
#define WINDOWRES_LEVELS 10
#define WINDOWRES_MAPS 100
#define WINDOWRES_PEAKS 100
#define WINDOWRES_PLOTS 10
#define WINDOWRES_PROFILES 100
#define WINDOWRES_COMPUTEMAX 10
#define WINDOWRES_PHASTCONS 10

/* Compact BG resolution                    */
#define COMPACTRES 10                    

/* FindPeaks default threshold              */
#define FINDPEAKS_T 1.5                    

/* Maximum gene length                      */
#define MAXGENELENGTH 2500000

/* Maximum size for promoter regions        */
#define MAXSIMPLIFIED 10000

/* Log10 noise threshold (heatmaps)         */
#define NOISETHRESHOLD 1

/* No ChIP signal                           */
#define NOSIGNAL 0

/* RefGene exon coordinates (L or R)        */
#define MAXEXONSINFO 5000
#define MAXEXONS 500

/* Pseudocount for logarithms               */
#define PSEUDOCOUNT 0.1

/* Extra space for READS bins               */
#define EXTRAMEM 100000

/* Million (for normalization)              */
#define MEGA 1000000

/* Input formats                            */
#define BEDGRAPH 0
#define WIG 1
#define SAM 2

/* Macros (functions)                       */
#define MIN(a,b) (a<b)?a:b;
#define MAX(a,b) (a>b)?a:b;


/*************************************************************************
B. DATA TYPES
*************************************************************************/

/* Information system */
typedef struct s_account               
{

  int nChrs;
  long nLines;

  long nTranscripts[MAXCHRS];
  long nTotalTranscripts;

  long nPeaks;
  long nPeaksA[MAXCHRS];
  long nPeaksB[MAXCHRS];
  long nTotalPeaksA;
  long nTotalPeaksB;

  long nReads;
  long nFwdReads;
  long nRvsReads;

  long nReads1;
  long nFwdReads1;
  long nRvsReads1;

  long nReads2;
  long nFwdReads2;
  long nRvsReads2;

  int nGenes;

  long nMemory;

  time_t tStart, tEnd;
} account;

/* Hash (dictionary) */
typedef struct s_node *pnode;          
typedef struct s_node
{
  char s[MAXSTRING];
  int key;
  pnode next;
} node; 

typedef struct s_dict                  
{
  pnode T[MAXENTRY];
  int nextFree; 
} dict;

/* Transcript information */
typedef struct s_transcript
{
  char chr[MAXSTRING];
  char strand;
  long pos1;
  long pos2;
  char gene[MAXGENENAME];
  char nm[MAXGENENAME];
} transcript; 

/* Peak information */
typedef struct s_peak
{
  char  chr[MAXSTRING];
  long  pos1;
  long  pos2;
  char  score[MAXSTRING];
  int   used;
} peak; 


/*************************************************************************
C. IMPORTED HEADERS
*************************************************************************/

/* Information system */
account* RequestMemoryAccounting();
account* InitAcc();

/* Output functions */
void printError(char *s);
void printMess(char* s);
void printRes(char* s);
void printSeqCodeHeader(char* s);
void StartingTime(char* mess, account* m);
void GetTime(char* mess, account* m);
void OutputTime(account *m);
void printReadingInit();
void printReadingInfo(long nLines);
void printReadingInfoSmall(long nLines);

void OutputMatch(char chr[],
		 long pos1A,
		 long pos2A,
		 char score[],
		 transcript t);

void OutputPeaksLists(peak** PeaksA,
		      peak** PeaksB,
		      char* OutputFileAB,
		      char* OutputFileA,
		      char* OutputFileB,
		      char* NameA,
		      char* NameB,
		      dict* ChrNames);

void OutputAvgReadsFile(long* REGION,
			int L,
			char* FileName,
			long ngenes,
			long nreads);

void OutputAvgReadsFileCombined(long* REGION,
				int L,
				char* FileName,
				long ngenes);

void OutputOneGeneHeatFile(char* gene,
                           long* REGION,
                           int L,
			   long nreads,
                           FILE* file);

void OutputOneGeneCombinedHeatFile(char* gene,
				   long* REGION,
				   int L,
				   FILE* file);

void OutputOneGeneGENEHeatFile(char* gene,
			       long* UPSTREAM,
			       long* IDEALGENE,
			       long* DOWNSTREAM,
			       int L,
			       long nreads,
			       FILE* file); 

void OutputIdealGeneProfileFile(long* UPSTREAM,
				long* IDEALGENE,
				long* DOWNSTREAM,
				int L,
				char* FileName,
				long nreads);

void OutputGeneAvgMax(transcript t, float avg, int max);

/* Hash functions */
void resetDict(dict* d);
int setkeyDict(dict* d, char s[]);
int getkeyDict(dict* d, char s[]);
void showDict(dict* d);

/* Read arguments */
void readargv_processmacs(int argc,
			  char* argv[], 
			  char* ChromInfoFile, 
			  char* MACSFile, 
			  char* OutputFile,
			  char* TrackName);

void readargv_matchpeaksgenes (int argc,
			       char* argv[],
			       char* PeaksFile, 
			       char* RefgeneFile,
			       int* Upstream,
			       int* Downstream,
			       int* TSSrange,
			       int* TESrange);

void readargv_matchpeaks(int argc,char* argv[],
			 char* PeaksFileA, 
			 char* PeaksFileB,
			 char* NameA,
			 char* NameB);

void readargv_produceTSSplots(int argc, char* argv[],
			      char ChromInfoFile[FILENAMELENGTH],
			      char RefgeneFile[FILENAMELENGTH],
			      char ProfileFile[FILENAMELENGTH],
			      char GenesFile[FILENAMELENGTH],
			      char TrackName[FILENAMELENGTH],
			      int* L);

void readargv_combineTSSmaps(int argc,char* argv[],
			     char ChromInfoFile[FILENAMELENGTH],
			     char RefgeneFile[FILENAMELENGTH],
			     char ProfileFile1[FILENAMELENGTH],
			     char ProfileFile2[FILENAMELENGTH],
			     char GenesFile[FILENAMELENGTH],
			     char TrackName1[FILENAMELENGTH],
			     char TrackName2[FILENAMELENGTH],
			     int* L);

void readargv_combineTSSplots(int argc,char* argv[],
			      char ChromInfoFile[FILENAMELENGTH],
			      char RefgeneFile[FILENAMELENGTH],
			      char ProfileFile1[FILENAMELENGTH],
			      char ProfileFile2[FILENAMELENGTH],
			      char GenesFile[FILENAMELENGTH],
			      char TrackName1[FILENAMELENGTH],
			      char TrackName2[FILENAMELENGTH],
			      int* L);

void readargv_produceGENEplots(int argc,char* argv[],
			       char ChromInfoFile[FILENAMELENGTH],
			       char RefgeneFile[FILENAMELENGTH],
			       char ProfileFile[FILENAMELENGTH],
			       char GenesFile[FILENAMELENGTH],
			       char TrackName[FILENAMELENGTH],
			       int* L);

void readargv_produceGENEmaps(int argc,char* argv[],
			      char ChromInfoFile[FILENAMELENGTH],
			      char RefgeneFile[FILENAMELENGTH],
			      char ProfileFile[FILENAMELENGTH],
			      char GenesFile[FILENAMELENGTH],
			      char TrackName[FILENAMELENGTH],
			      int* L);

void readargv_computemaxsignal(int argc,char* argv[],
			       char ChromInfoFile[FILENAMELENGTH],
			       char RefgeneFile[FILENAMELENGTH],
			       char ProfileFile[FILENAMELENGTH],
			       int* Upstream);

void readargv_buildChIPprofile(int argc,char* argv[],
			       char ChromInfoFile[FILENAMELENGTH],
			       char ProfileFile[FILENAMELENGTH],
			       char TrackName[FILENAMELENGTH]);

void readargv_combineChIPprofiles(int argc,char* argv[],
				   char ChromInfoFile[FILENAMELENGTH],
				   char ProfileFile1[FILENAMELENGTH],
				   char ProfileFile2[FILENAMELENGTH],
				   char TrackName1[FILENAMELENGTH],
				   char TrackName2[FILENAMELENGTH]);

void readargv_findPeaks(int argc,char* argv[],
			char ChromInfoFile[FILENAMELENGTH],
			char ProfileFile[FILENAMELENGTH],
			char TrackName[FILENAMELENGTH]);

void readargv_genomeDistribution(int argc,char* argv[],
				 char ChromInfoFile[FILENAMELENGTH],
				 char RefgeneFile[FILENAMELENGTH],
				 char Peaks[FILENAMELENGTH],
				 char TrackName[FILENAMELENGTH]);

void readargv_produceTESplots(int argc, char* argv[],
			      char ChromInfoFile[FILENAMELENGTH],
			      char RefgeneFile[FILENAMELENGTH],
			      char ProfileFile[FILENAMELENGTH],
			      char GenesFile[FILENAMELENGTH],
			      char TrackName[FILENAMELENGTH],
			      int* L);

void readargv_producePEAKplots(int argc,char* argv[],
			       char ChromInfoFile[FILENAMELENGTH],
			       char ReadsFile[FILENAMELENGTH],
			       char PeaksFile[FILENAMELENGTH],
			       char TrackName[FILENAMELENGTH],
			       int* L);

void readargv_recoverChIPlevels(int argc,char* argv[],
				char ChromInfoFile[FILENAMELENGTH],
				char ReadsFile[FILENAMELENGTH],
				char PeaksFile[FILENAMELENGTH],
				char TrackName[FILENAMELENGTH]);

void readargv_scorePhastCons(int argc,char* argv[],
			     char ChromInfoFile[FILENAMELENGTH],
			     char RefgeneFile[FILENAMELENGTH],
			     char ProfileFile[FILENAMELENGTH],
			     char PeaksFile[FILENAMELENGTH],
			     char TrackName[FILENAMELENGTH],
			     float* Threshold);

void readargv_produceTSSmaps(int argc,char* argv[],
			     char ChromInfoFile[FILENAMELENGTH],
			     char RefgeneFile[FILENAMELENGTH],
			     char ProfileFile[FILENAMELENGTH],
			     char GenesFile[FILENAMELENGTH],
			     char TrackName[FILENAMELENGTH],
			     int* L);

void readargv_produceTESmaps(int argc,char* argv[],
			     char ChromInfoFile[FILENAMELENGTH],
			     char RefgeneFile[FILENAMELENGTH],
			     char ProfileFile[FILENAMELENGTH],
			     char GenesFile[FILENAMELENGTH],
			     char TrackName[FILENAMELENGTH],
			     int* L);

void readargv_producePEAKmaps(int argc,char* argv[],
			      char ChromInfoFile[FILENAMELENGTH],
			      char ReadsFile[FILENAMELENGTH],
			      char PeaksFile[FILENAMELENGTH],
			      char TrackName[FILENAMELENGTH],
			      int* L);

/* Request memory */
dict* RequestMemoryDictionary();
long* RequestMemoryChrSizes();
transcript** RequestMemoryTranscripts();
peak** RequestMemoryPeaks();
int* RequestMemoryGeneChrs();
unsigned int** RequestMemoryREADS(long* ChrSizes, dict* ChrNames);
long* RequestMemoryRegion(int L);
float** RequestMemoryPHASTCONS(long* ChrSizes, dict* ChrNames);

/* Reading input files */
long GetFileSize(char* File);

int ReadChrFile (char *FileName, long* ChrSizes, dict* ChrNames);

long ReadRefgeneFile(char* RefgeneFile,
		     transcript** Transcripts,
		     int Upstream,
		     int Downstream,
		     int TSSrange,
		     int TESrange,
		     dict* ChrNames);

long ReadPeaksFile(char* FileName,
		   peak** Peaks, 
		   long* nPeaksX,
		   dict* ChrNames);

long ReadRefgeneFile2(char* FileName,
		      transcript** Transcripts,
		      dict* ChrNames,
		      dict* GeneNames,
		      int* GeneChrs);

void ReadSAMFile (char *FileName, 
		  long* ChrSizes, 
		  dict* ChrNames, 
		  unsigned int** READS);

void ReadSAMPairedFile (char *FileName, 
			long* ChrSizes, 
			dict* ChrNames, 
			unsigned int** READS);

long ReadBGFile (char *FileName, 
                 long* ChrSizes, 
                 dict* ChrNames, 
                 unsigned int** READS);

long ReadRefgeneFileMark (char *FileName, 
			  long* ChrSizes, 
			  dict* ChrNames, 
			  unsigned int** REGIONS);

long ReadRefgeneFileMarkExons (char *FileName, 
			       long* ChrSizes, 
			       dict* ChrNames, 
			       unsigned int** EXONS);

void ReadPhastConsFiles (char* Folder, 
			 long* ChrSizes, 
			 dict* ChrNames, 
			 float** PHASTCONS);

long ReadWIGFile (char *FileName, 
		  long* ChrSizes, 
		  dict* ChrNames, 
		  float** PHASTCONS);

/* Processing routines */
long ReadProcessPeaks(char* FileName,
                      long* ChrSizes,
                      dict* ChrNames,
                      int L,
                      unsigned int** READS,
                      long* REGION,
                      char* OutputFileName,
                      char* OutputFileName2,
                      long nreads);

long ProcessBGFile (char* FileName, 
		    long* ChrSizes, 
		    dict* ChrNames,
		    char* OutputFile,
		    char* TrackName);

long ProcessPeaks(char* PeaksFile,transcript** Transcripts, dict* ChrNames);

void ProcessPeaksLists(peak** PeaksA, 
		       peak** PeaksB, 
		       dict* ChrNames,
		       char* NameA);

int ProcessTargetGenes(char* FileName,
                       transcript** Transcripts,
                       dict* ChrNames,
		       long* ChrSizes,
                       dict* GeneNames,
                       int* GeneChrs,
                       int L,
                       unsigned int** READS,
                       long * REGION,
		       char* OutputFileName,
		       char* OutputFileName2,
                       long nreads);

int CombineTargetGenes(char* FileName,
		       transcript** Transcripts, 
		       dict* ChrNames,
		       long* ChrSizes,
		       dict* GeneNames,
		       int* GeneChrs,
		       int L,
		       unsigned int** READS1,
		       unsigned int** READS2,
		       long* REGION,
		       char* OutputFileName,
		       long nreads1,
		       long nreads2);

int CombineTargetGenesMaps(char* FileName,
			   transcript** Transcripts, 
			   dict* ChrNames,
			   long* ChrSizes,
			   dict* GeneNames,
			   int* GeneChrs,
			   int L,
			   unsigned int** READS1,
			   unsigned int** READS2,
			   char* OutputFileName2,
			   long nreads1,
			   long nreads2);

int ProcessTargetGenesGENE(char* FileName,
			   transcript** Transcripts, 
			   dict* ChrNames,
			   long* ChrSizes,
			   dict* GeneNames,
			   int* GeneChrs,
			   int L,
			   unsigned int** READS,
			   long* UPSTREAM,
			   long* IDEALGENE,
			   long* DOWNSTREAM,
			   char* OutputFileName2,
			   long nreads);

int ProcessGenesMax(char RefgeneFile[FILENAMELENGTH],
		    transcript** Transcripts, 
		    dict* ChrNames,
		    long* ChrSizes,
		    dict* GeneNames,
		    int* GeneChrs,
		    int Upstream,
		    unsigned int** READS);

int ProcessGenome(dict* ChrNames,
		  long* ChrSizes,
		  unsigned int** READS,
		  char* TrackName,
		  char* TrackColor,
		  char* OutputFileName,
		  long nreads);

int CombineGenomes(dict* ChrNames,
		   long* ChrSizes,
		   unsigned int** READS1,
		   unsigned int** READS2,
		   char* TrackName1,
		   char* TrackName2,
		   char* TrackColor,
		   char* OutputFileName,
		   long nreads1,
		   long nreads2);

int ProcessGenomePeaks(dict* ChrNames,
		       long* ChrSizes,
		       unsigned int** READS,
		       char* TrackName,
		       char* OutputFileName,
		       long nreads);

long ProcessPeaksMarks(char* FileName,
		       dict* ChrNames,
		       long* ChrSizes,  
		       unsigned int** REGIONS,
		       char* OutputFileName,
		       char* RscriptFileName,
		       char* PlotFileName,
		       char* TrackName);

int ProcessTargetGenesTES(char* FileName,
			  transcript** Transcripts,
			  dict* ChrNames,
			  long* ChrSizes,
			  dict* GeneNames,
			  int* GeneChrs,
			  int L,
			  unsigned int** READS,
			  long * REGION,
			  char* OutputFileName,
			  long nreads);

long ReadProcessPeaksLevels(char* FileName,
			    long* ChrSizes,
			    dict* ChrNames,
			    unsigned int** READS,
			    char* OutputFileName,
			    long nreads);

long ReadProcessPeaksPhastCons(char* FileName,
			       long* ChrSizes,
			       dict* ChrNames,
			       float** PHASTCONS,
			       unsigned int** EXONS,
			       char* OutputFileName);

int ProcessTargetGenesMaps(char* FileName,
                           transcript** Transcripts,
                           dict* ChrNames,
                           long* ChrSizes,
                           dict* GeneNames,
                           int* GeneChrs,
                           int L,
                           unsigned int** READS,
                           char* OutputFileName2,
                           long nreads);

int ProcessTargetGenesGENEMaps(char* FileName,
			       transcript** Transcripts, 
			       dict* ChrNames,
			       long* ChrSizes,
			       dict* GeneNames,
			       int* GeneChrs,
			       int L,
			       unsigned int** READS,
			       char* OutputFileName2,
			       long nreads);

int ProcessTargetGenesTESMaps(char* FileName,
			      transcript** Transcripts,
			      dict* ChrNames,
			      long* ChrSizes,
			      dict* GeneNames,
			      int* GeneChrs,
			      int L,
			      unsigned int** READS,
			      char* OutputFileName2,
			      long nreads);

long ReadProcessPeaksMaps(char* FileName,
			  long* ChrSizes,
			  dict* ChrNames,
			  int L,
			  unsigned int** READS,
			  char* OutputFileName2,
			  char* OutputFileName3,
			  long nreads);

/* R script generation */
void CreateRscript(char* OutputFileName,
		   char* RscriptFileName,
		   char* PlotFileName,
		   char* TrackName,
		   int L,
		   long nreads,
		   long nGenes);

void CreateRscriptCombined(char* OutputFileName,
                           char* RscriptFileName,
                           char* PlotFileName,
                           char* TrackName1,
                           char* TrackName2,
                           int L,
                           long nreads1,
                           long nreads2,
                           long nGenes);

void CreateRscriptCombinedMaps(char* OutputFileName2,
			       char* RscriptFileName,
			       char* PlotFileName2,
			       char* TrackName1,
			       char* TrackName2,
			       int L,
			       long nreads1,
			       long nreads2,
			       long nGenes);

void CreateRscriptGENE(char* OutputFileName,
		       char* RscriptFileName,
		       char* PlotFileName,
		       char* TrackName,
		       int L,
		       long nreads);

void CreateRscriptChart(char* RscriptFileName,
			char* PlotFileName,
			int nPromoters,
			int nIntragenics,
			int nIntergenics,
			char* TrackName,
			long nPeaks);

void CreateRscriptChartSimplified(char* RscriptFileName,
				  char* PlotFileName,
				  int nPromoters,
				  int nIntragenics,
				  int nIntergenics,
				  char* TrackName,
				  long nPeaks);

void CreateRscriptChartDetailed(char* RscriptFileName,
				char* PlotFileName,
				int nUtr5s,
				int nUtr3s,
				int nCDSs,
				int nDistals,
				int nProximals,
				int nIntrons,
				int nIntergenics,
				char* TrackName,
				long nPeaks);

void CreateRscriptTES(char* OutputFileName,
		      char* RscriptFileName,
		      char* PlotFileName,
		      char* TrackName,
		      int L,
		      long nreads,
		      long nGenes);

void CreateRscriptPeak(char* OutputFileName,
		       char* RscriptFileName,
		       char* PlotFileName,
		       char* TrackName,
		       int L,
		       long nreads,
		       long nPeaks);

void CreateRscriptMaps(char* OutputFileName2,
		       char* RscriptFileName,
		       char* PlotFileName2,
		       char* TrackName,
		       int L,
		       long nreads,
		       long nGenes);

void CreateRscriptGENEMaps(char* OutputFileName2,
			   char* RscriptFileName,
			   char* PlotFileName2,
			   char* TrackName,
			   int L,
			   long nreads,
			   long nGenes);

void CreateRscriptTESMaps(char* OutputFileName2,
			  char* RscriptFileName,
			  char* PlotFileName2,
			  char* TrackName,
			  int L,
			  long nreads,
			  long nGenes);


void CreateRscriptPeakMaps(char* OutputFileName4,
			   char* RscriptFileName,
			   char* PlotFileName2,
			   char* TrackName,
			   int L,
			   long nreads,
			   long nPeaks);

/* R color selection */

void LoadRColors(dict* RColors,
		 char RGBcodes[MAXCOLORS][MAXRGBCODE]);

void ValidateRColors(char BGCOLOR[MAXCOLORNAME],
		     char FGCOLOR[MAXCOLORNAME],
		     char GBGCOLOR[MAXCOLORNAME],
		     char GFGCOLOR[MAXCOLORNAME],
		     char HBGCOLOR[MAXCOLORNAME],
		     char FGCOLORB[MAXCOLORNAME],
		     dict* RColors);

void ValidateOneRColor(char FGCOLOR[MAXCOLORNAME],
		       dict* RColors,
		       char RGBcodes[MAXCOLORS][MAXRGBCODE],
		       char* TrackColor);

void ValidateRColorsPieCharts(char COLOR1[MAXCOLORNAME],
			      char COLOR2[MAXCOLORNAME],
			      char COLOR3[MAXCOLORNAME],
			      char COLOR4[MAXCOLORNAME],
			      char COLOR5[MAXCOLORNAME],
			      char COLOR6[MAXCOLORNAME],
			      char COLOR7[MAXCOLORNAME],
			      dict* RColors);

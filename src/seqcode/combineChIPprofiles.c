/*************************************************************************
*                                                                        *
*   Module: combineChIPprofiles                                          *
*                                                                        *
*   CombineChIPprofiles main program                                     *
*                                                                        *
*   This file is part of the SeqCode 1.0 distribution                    *
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

#include "seqcode/seqcode.h"

/************************************************************************
                            ENVIRONMENTAL variables
************************************************************************/

/* Verbose flag (memory/processing information) */
int VRB=0;

/* Read extension */
int AVGL = 150;

/* Bin resolution */
int WINDOWRES = WINDOWRES_PROFILES;

/* Spike-in normalization */
int SPIKEIN = 0;

/* Pair-end processing */
int PAIREND = FALSE;

/* Minus-strand (strand-specific RNAseq scores) */
int MINUSSTRAND = FALSE;

/* Track color */
char FGCOLOR[MAXCOLORNAME];

/* Accounting time and results */
account *m;

/************************************************************************
                            MAIN program
************************************************************************/

int main (int argc, char *argv[])
{
  char ChromInfoFile[FILENAMELENGTH],
    ReadsFile1[FILENAMELENGTH],
    ReadsFile2[FILENAMELENGTH],
    TrackName1[FILENAMELENGTH],
    TrackName2[FILENAMELENGTH];
  
  long* ChrSizes;
  dict* ChrNames;

  unsigned int** READS1;
  unsigned int** READS2;

  long FileSize;

  char OutputFileName[FILENAMELENGTH];
  char OutputDirName[FILENAMELENGTH];

  int nChromosomes;

  char tmpFileName[FILENAMELENGTH];

  dict* RColors;
  char RGBcodes[MAXCOLORS][MAXRGBCODE];

  char TrackColor[MAXRGBCODE];

  int tmpSystem;
  
  char mess[MAXSTRING];

  
  /** 0. Starting and reading options, parameters and sequence... **/
  
  /* 0.a. Initializing stats and time counters */
  m = (account *) InitAcc();  
  
  /* 0.b. Read setup options */
  readargv_combineChIPprofiles(argc,argv,
			    ChromInfoFile,
			    ReadsFile1,
			    ReadsFile2,
			    TrackName1,
			    TrackName2);
  printSeqCodeHeader(COMBINECHIPPROFILES);
  StartingTime(mess,m);
  printMess(mess);

  /* 1. Allocating data structures */
  printMess("1. Request Memory to Operating System");
  ChrNames = (dict *) RequestMemoryDictionary();
  ChrSizes = (long *) RequestMemoryChrSizes();
  RColors = (dict *) RequestMemoryDictionary();

  /* Color control */
  LoadRColors(RColors,RGBcodes);
  ValidateOneRColor(FGCOLOR,RColors,RGBcodes,TrackColor);
  printRes("Loading the R color schema");

  /* 2. Read the ChrSize file */
  printMess("2. Reading Chromosome Sizes");
  m->nChrs = ReadChrFile (ChromInfoFile,ChrSizes,ChrNames);
  sprintf(mess,"Size was successfully acquired for %d chromosomes",m->nChrs);
  printRes(mess);

  /*  Additional allocating data structures */
  READS1 = (unsigned int**) RequestMemoryREADS(ChrSizes,ChrNames);
  READS2 = (unsigned int**) RequestMemoryREADS(ChrSizes,ChrNames);

  /* 3. Create the output files (BedGraph ChIPseq/RNAseq profile */
  printMess("3. Creating output file");
  
  sprintf(OutputDirName,"%s-%s_BedGraph_profile",TrackName1,TrackName2);
  sprintf(mess,"Creating the output folder %s",OutputDirName);
  printRes(mess);

  /* Creating the new Output folder */
  if (mkdir(OutputDirName,0777))
    {
      printRes("Removing the old output folder");
      sprintf(mess,"rm -rf %s",OutputDirName);
      printRes(mess);
      tmpSystem = system(mess);
      
      if (mkdir(OutputDirName,0777))
	{
	  printError("Problems with the creation of the output folder");
	}
    }

  /* BedGraph file */
  sprintf(OutputFileName,"%s/%s-%s.bedgraph",
	  OutputDirName,
	  TrackName1,
	  TrackName2);
  sprintf(mess,"Filename for the BedGraph profile %s",OutputFileName);
  printRes(mess);

  /* 4a. Acquire the BAM/SAM reads file 1 */
  printMess("4a. Processing the BAM/SAM file 1");
  
  FileSize = GetFileSize(ReadsFile1);
  sprintf(mess,"%s: %.2lf Mb",ReadsFile1,(float)FileSize/MEGABYTE);
  printRes(mess);

  if (!PAIREND)
    {
      printRes("Processing Single-end sequencing reads");

      ReadSAMFile(ReadsFile1,ChrSizes,ChrNames,READS1);
      m->nReads1 = m->nReads;
      m->nFwdReads1 = m->nFwdReads;
      m->nRvsReads1 = m->nRvsReads;
      sprintf(mess,"%ld reads of the SAM file processed",m->nReads1);
      printRes(mess);
      sprintf(mess,"%ld forward reads acquired",m->nFwdReads1);
      printRes(mess);
      sprintf(mess,"%ld reverse reads acquired",m->nRvsReads1);
      printRes(mess);
    }
  else
    {
      printRes("Processing Pair-end sequencing reads");

      ReadSAMPairedFile(ReadsFile1,ChrSizes,ChrNames,READS1);
      m->nReads1 = m->nReads;
      m->nFwdReads1 = m->nFwdReads;
      m->nRvsReads1 = m->nRvsReads;
      sprintf(mess,"%ld reads of the SAM file processed",m->nReads1);
      printRes(mess);
      sprintf(mess,"%ld reads mapped in the proper pair",
	      m->nFwdReads1 + m->nRvsReads1);
      printRes(mess);
      sprintf(mess,"%ld forward reads acquired",m->nFwdReads1);
      printRes(mess);
      sprintf(mess,"%ld reverse reads acquired",m->nRvsReads1);
      printRes(mess);

    }

  /* 4b. Acquire the BAM/SAM reads file 2 */
  printMess("4b. Processing the BAM/SAM file 2");
  
  FileSize = GetFileSize(ReadsFile2);
  sprintf(mess,"%s: %.2lf Mb",ReadsFile2,(float)FileSize/MEGABYTE);
  printRes(mess);

  if (!PAIREND)
    {
      printRes("Processing Single-end sequencing reads");

      ReadSAMFile(ReadsFile2,ChrSizes,ChrNames,READS2);
      m->nReads2 = m->nReads;
      m->nFwdReads2 = m->nFwdReads;
      m->nRvsReads2 = m->nRvsReads;
      sprintf(mess,"%ld reads of the SAM file processed",m->nReads2);
      printRes(mess);
      sprintf(mess,"%ld forward reads acquired",m->nFwdReads2);
      printRes(mess);
      sprintf(mess,"%ld reverse reads acquired",m->nRvsReads2);
      printRes(mess);
    }
  else
    {
      printRes("Processing Pair-end sequencing reads");

      ReadSAMPairedFile(ReadsFile2,ChrSizes,ChrNames,READS2);
      m->nReads2 = m->nReads;
      m->nFwdReads2 = m->nFwdReads;
      m->nRvsReads2 = m->nRvsReads;
      sprintf(mess,"%ld reads of the SAM file processed",m->nReads2);
      printRes(mess);
      sprintf(mess,"%ld reads mapped in the proper pair",
	      m->nFwdReads2+m->nRvsReads2);
      printRes(mess);
      sprintf(mess,"%ld forward reads acquired",m->nFwdReads2);
      printRes(mess);
      sprintf(mess,"%ld reverse reads acquired",m->nRvsReads2);
      printRes(mess);

    }

  /* 5. Process the chromosomes to combine the BedGraph profile */
  printMess("5. Combining the BedGraph profile");

  nChromosomes = CombineGenomes(ChrNames,
				ChrSizes,
				READS1,
				READS2,
				TrackName1,
				TrackName2,
				TrackColor,
				OutputFileName,
				m->nReads1,
				m->nReads2);
  sprintf(mess,"%d chromosomes sucessfully included into the BG profile\n",
	  nChromosomes);
  printRes(mess);

  /* Inform about the final output file size (before zip) */
  printRes("Size of the resulting BG file (uncompressed):");
  FileSize = GetFileSize(OutputFileName);
  sprintf(mess,"%s: %.2lf Mb",OutputFileName,(float)FileSize/MEGABYTE);
  printRes(mess);

  /* Compress the resulting BedGraph file */
  printRes("Compressing the final BG profile\n");
  sprintf(mess,"gzip %s",OutputFileName);
  printRes(mess);
  tmpSystem = system(mess);
  
  /* Inform about the final output file size (after zip) */
  printRes("Size of the resulting BG file:");
  sprintf(tmpFileName,"%s.gz",OutputFileName);
  FileSize = GetFileSize(tmpFileName);
  sprintf(mess,"%s: %.2lf Mb",tmpFileName,(float)FileSize/MEGABYTE);
  printRes(mess);

  /* 6. The End */
  OutputTime(m);  
  
  exit(0);
  return(0);
}

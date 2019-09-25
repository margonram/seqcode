/*************************************************************************
*                                                                        *
*   Module: readargv_combineChIprofiles                                  *
*                                                                        *
*   Read set up options and filenames from user input (buildChIPprofile) *
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

/* SeqCode external vars */
extern int VRB;
extern int AVGL;
extern int WINDOWRES;
extern int PAIREND;
extern int MINUSSTRAND;
extern char FGCOLOR[MAXCOLORNAME];

/* required by getopts */
extern char* optarg;
extern int optind;

/* usage string */
char USAGE[MAXUSAGE][MAXSTRING];

void printHelp_combineChIPprofiles()
{
  fprintf(stderr,"%s\n",USAGE[0]);
  fprintf(stderr,"%s\n",USAGE[1]);
  fprintf(stderr,"%s\n",USAGE[2]);
  fprintf(stderr,"%s\n",USAGE[3]);

  fprintf(stderr,"OPTIONS\n");
  fprintf(stderr,"\t-c : Track color (default: black).\n");
  fprintf(stderr,"\t-l : Avg. fragment size (default: 150).\n");
  fprintf(stderr,"\t-m : Invert the score of BedGraph records.\n");
  fprintf(stderr,"\t-p : Using pair-end reads mapped in proper pair (default: single-end).\n");
  fprintf(stderr,"\t-w : Window resolution (default: 100).\n");
  fprintf(stderr,"\t-v : Verbose. Display info messages.\n");
  fprintf(stderr,"\t-h : Show this help.\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"%s\n",USAGE[4]);
  fprintf(stderr,"%s\n",USAGE[5]);
  fprintf(stderr,"%s\n",USAGE[6]);
}

void readargv_combineChIPprofiles (int argc,char* argv[],
				   char ChromInfoFile[FILENAMELENGTH],
				   char ProfileFile1[FILENAMELENGTH],
				   char ProfileFile2[FILENAMELENGTH],
				   char TrackName1[FILENAMELENGTH],
				   char TrackName2[FILENAMELENGTH])
{
  int c;
  char mess[MAXSTRING];

  int nFiles;

  /* Adding usage messages */
  sprintf(USAGE[0],
	  "%s\t\t\t\t\t\tUser commands\t\t\t\t\t\t%s\n",
	  VERSION,
	  COMBINECHIPPROFILES);
  sprintf(USAGE[1],
	  "NAME\n\t%s - a program to combine two genome-wide ChIPseq/ATACseq/RNAseq profiles.\n",
	  COMBINECHIPPROFILES);
  sprintf(USAGE[2],
	  "SYNOPSIS\n\t%s [-c <Rcolor>][-l <bp>][-m][-p][-w <bp>][-v][-h]\n\t<chrom_info> <SAM file1> <SAM file2> <name1> <name2>\n",
	  COMBINECHIPPROFILES);
  sprintf(USAGE[3],
	  "OUTPUT\n\tOne folder with one file:\n\t- BedGraph profile (compressed).\n");
  sprintf(USAGE[4],
	  "SEE ALSO\n\tSeqCode homepage: %s\n\tGitHub source code: %s\n",
	  HOMEPAGE,SOURCECODE);
  sprintf(USAGE[5],
	  "AUTHORS\n\tWritten by %s.\n",
	  AUTHORSHIP);
  sprintf(USAGE[6],
	  "%s\t\t\t\t\t\tUser commands\t\t\t\t\t\t%s\n",
	  VERSION,
	  COMBINECHIPPROFILES);

  /* Default background/foreground colors */
  strcpy(FGCOLOR,"black");
  
  /* Reading setup options */
  while ((c = getopt(argc,argv,"c:l:mpw:vh")) != -1)
    switch(c)
      {
      case 'c':
        strcpy(FGCOLOR,optarg);
        break;
      case 'l':
	AVGL = atoi(optarg);
	break;
      case 'm': MINUSSTRAND++;
	break;
      case 'p': PAIREND++;
	break;
      case 'w':
        WINDOWRES = atoi(optarg);
        break;
      case 'v': VRB++;
	break;
      case 'h': printHelp_combineChIPprofiles();
	exit(1);
	break;
      }
  
  /* Setup Errors : Wrong number of filenames */
  nFiles = argc - optind;

  if (nFiles != 5)
    {
      fprintf(stderr,"%s\n",USAGE[0]);
      fprintf(stderr,"%s\n",USAGE[1]);
      fprintf(stderr,"%s\n",USAGE[2]);

      sprintf(mess,"Five parameters are expected, %d are provided instead\n",
	      nFiles);
      printError(mess);
    }

  /* Acquire the name of the input files  */ 
  strcpy(ChromInfoFile,argv[optind]); 
  strcpy(ProfileFile1,argv[optind+1]);
  strcpy(ProfileFile2,argv[optind+2]);
  strcpy(TrackName1,argv[optind+3]);
  strcpy(TrackName2,argv[optind+4]);
}

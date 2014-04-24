/** Testing of read mapping strategies in rmap.c
 */
/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                           * 
 *                                                                           *        
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                              *
 *                                                                           *
 *  This file is part of SMALT.                                              *
 *                                                                           *
 *  SMALT is free software: you can redistribute it and/or modify it under   *
 *  the terms of the GNU General Public License as published by the Free     *
 *  Software Foundation, either version 3 of the License, or (at your        *
 *  option) any later version.                                               *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU         *
 *  General Public License for more details.                                 *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

#include "elib.h"
#include "sequence.h"
#include "hashidx.h"
#include "results.h"
#include "resultpairs.h"
#include "report.h"
#include "rmap.h"

enum {
  RMAP_MIN_BASQUAL = 0,
};

typedef unsigned char UCHAR;
typedef unsigned int UINT32;

static int calcMinKtup(int mincov_percent, int readlen, UCHAR ktup, UCHAR nskip)
{
  int n_tup = (readlen-ktup+1)/nskip;
  int minktup = (n_tup-1)*mincov_percent/100;
  return ((minktup))? minktup: 1;
}

int main(int argc, char *argv[])
{
  int errcode = ERRCODE_SUCCESS;
  char is_verbose = 0;
  char *binfilnam, *readfilnam;
  UCHAR nskip, rmapflg = 0, rsltflg = 0;
  REPMODIFLG_t repmodflg = 0;
  REPOUFMT_t repfmt = REPORTFMT_CIGAR;
  short max_depth = 4000, target_depth;
  int maxhit_per_tuple, depthtarget, mincov_percent;
  int min_swatscor = 13;
  int minswatscor_below_max = 20;
  int sctr;
  UINT32 readlen, min_cover;
  SeqCodec *codecp;
  SeqFastq *readp, *rbufp, *qbufp;
  SeqSet *ssp;
  SeqIO *sfp;
  HashTable *htp;
  ScorePenalties *scorpltyp;
  ScoreMatrix *scormtxp;
  ResultFilter *rfp = resultSetCreateFilter();
  const ResultSet *rsltp;
  Report *rep;
  ReportWriter *repwp;
  RMap *rmp;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc != 8) {
    printf("usage: %s <hash table root name> <ncut> <target_depth> <mincover [%%]>" \
	   " <best only [y/n]> <alignment output [y/n]> <FASTQ file>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  
  binfilnam = argv[1];
  maxhit_per_tuple = atoi(argv[2]);
  depthtarget = atoi(argv[3]);
  if (depthtarget > SHRT_MAX) {
    printf("<target_depth> %i has to be a number <= %hi\n", depthtarget, SHRT_MAX);
    exit(EXIT_FAILURE);
  }
  target_depth = (short) depthtarget;
    
  mincov_percent = atoi(argv[4]);
  if (toupper(argv[5][0]) == 'Y' &&  argv[5][1] == '\0') {
    rsltflg |= RESULTFLG_BEST;
    minswatscor_below_max = 0;
  }
  if (toupper(argv[6][0]) == 'Y' &&  argv[6][1] == '\0')
    repmodflg |= REPORTMODIF_ALIOUT;
  readfilnam = argv[7];

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (NULL == (scorpltyp = scorePenaltiesCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (NULL == (scormtxp = scoreCreateMatrix(codecp, scorpltyp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(readp = seqFastqCreate(0, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  resultSetFilterData(rfp, min_swatscor, minswatscor_below_max, .0);

  /* load data */
  printf("Reading reference sequences ...\n");
  ssp = seqSetReadBinFil(&errcode, binfilnam);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);

  printf("Reading hash table ...\n");
  htp = hashTableRead(&errcode, binfilnam);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);
  hashTablePrintStats(stdout, htp);
  hashTableGetKtupLen(htp, &nskip);
  sfp = seqIOopen(&errcode, readfilnam, SEQIO_READ, 0);
  if (errcode) ERRMSGNO(errmsgp, errcode);

  if (!(rmp = rmapCreate(htp, codecp, ssp, NULL, 0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(rep = reportCreate(0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(repwp = reportCreateWriter(&errcode, "-", repfmt, repmodflg, ssp,
				   "rmap_test", NULL, argv, argc)))
    ERRMSGNO(errmsgp, errcode);

  printf("Working query sequence files ...\n");
  sctr = 0;
  while(!seqIOstatus(sfp)) {
    if ((is_verbose))
      printf("Reading query sequence %i ...\n", sctr++);
    if ((errcode = seqFastqRead(readp, sfp)))
      ERRMSGNO(errmsgp, errcode);
    if ((errcode = seqFastqEncode(readp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    seqFastqGetConstSequence(readp, &readlen, NULL);
    if (readlen > INT_MAX)
      ERRMSG(errmsgp, "query read", ERRCODE_LONGSEQ);
    min_cover = readlen*mincov_percent/100;

    rmapSingle(errmsgp,
	       rmp, 
#ifdef hashhit_dump_sortarray
	       NULL,
	       NULL,
#endif
	       readp, maxhit_per_tuple,
	       min_cover, min_swatscor, minswatscor_below_max, 
	       RMAP_MIN_BASQUAL,
	       target_depth, max_depth,
	       rmapflg, scormtxp,
	       rfp, htp, ssp, codecp);
    rmapGetData(&rsltp, NULL, NULL, &rbufp, &qbufp, rmp);
    errcode = resultSetAddToReport(rep, rsltflg, rsltp);
        
    if (errcode) 
      ERRMSGNO(errmsgp, errcode);
      errcode = reportWrite(repwp, 
			    readp, NULL,
			    ssp, codecp,
			    rep);
  }
  if (seqIOstatus(sfp) && 
      seqIOstatus(sfp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp));
  seqIOclose(sfp);

  /* clean up */
  hashTableDelete(htp);
  seqSetDelete(ssp);
  reportDeleteWriter(repwp);
  reportDelete(rep);
  rmapDelete(rmp);
  seqFastqDelete(readp);
  scoreDeleteMatrix(scormtxp);
  scorePenaltiesDelete(scorpltyp);
  seqCodecDelete(codecp);
  resultSetDeleteFilter(rfp);
  ERRMSG_END(errmsgp);

  return ERRCODE_SUCCESS;
}

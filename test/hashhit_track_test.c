/** Test seed lookup from hash table and check
 * coverage of simulated reads.
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

#include "elib.h"
#include "sequence.h"
#include "hashidx.h"
#include "hashhit.h"

enum {
  HASH_SEQREAD_BUFSIZ = 128*1024*1024,
  SEQBLOCKSIZ = 64*1024*1024,
  MAX_NHIT_PER_TUPLE = 10000,
  MAX_NHIT_TOTAL = 32000,
  HASH_MAXN_HITS = 1024*1024,
  HITINFO_BLKSZ =  1024,
  HASH_MIN_BASQUAL = 0,
  TRACKER_IS1BASED = 1,
};

typedef unsigned char UCHAR;

#ifdef hashhit_dump_sortarray
static const char hashhit_helper_output_filnam[] = "hashhit_track.bin";
#endif

static int processRead(ErrMsg *errmsgp, 
		       UCHAR is_verbose,
#ifdef RESULTS_TRACKER
		       Track *trkp,
#endif
#ifdef hashhit_dump_sortarray
		       FILE *dumpfp,
		       int *dumpctr,
#endif
		       HashHitInfo *hhip, HashHitList *hhlp,
		       const HashTable *htp, const SeqSet *ssp,
		       const SeqFastq *sqp, const SeqCodec *codecp)
{
  int errcode;
  unsigned char is_reverse;
  const SETSIZ_t *soffsp;
  const SEQNUM_t nseq = seqSetGetOffsets(ssp, &soffsp);

#ifdef RESULTS_TRACKER
  if ((errcode = trackMakeFromSequence(trkp, sqp, TRACKER_IS1BASED, ssp, htp)))
     ERRMSGNO(errmsgp, errcode);
#endif

  for (is_reverse = 0; is_reverse <= 1; is_reverse++) {
    SEQNUM_t s;
    if (is_verbose)
      printf("Get hit info %s...\n", (is_reverse)? "reverse":"forward");

    if ((errcode = hashCollectHitInfoShort(hhip, is_reverse, MAX_NHIT_PER_TUPLE, 
					   MAX_NHIT_TOTAL, HASH_MIN_BASQUAL, 
					   sqp, htp)))
      ERRMSGNO(errmsgp, errcode);

    
    for (s=0; s<nseq; s++) {
      hashBlankHitList(hhlp);
      if ((errcode = hashCollectHitsForSegment(hhlp, 
#ifdef RESULTS_TRACKER
					       trkp,
#endif
#ifdef hashhit_dump_sortarray
					       dumpfp,
					       dumpctr,
#endif
					       soffsp[s], soffsp[s+1],
#ifndef hashhit_minimise_coverdeficit 
					       MAX_NHIT_PER_TUPLE,
					       1,
#endif
					       hhip, htp, NULL)))
	break;
      if (is_verbose)
	printf("Checking list ...\n");

      if ((errcode = hashCheckHitList(hhlp, sqp, htp, ssp, codecp)))
	ERRMSGNO(errmsgp, errcode);
    }
  }
#ifdef RESULTS_TRACKER
  if (is_verbose)
    trackPrintf(stdout, trkp);
#endif
  return errcode;
}

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  int sctr;
  char *binfilnam, *readfilnam;
  UCHAR is_verbose = 0;
  SeqIO *sfp;
  SeqCodec *codecp;
  SeqFastq *sqp;
  SeqSet *ssp;
  HashTable *htp;
  HashHitInfo *hhip;
  HashHitList *hhlp;
  ErrMsg *errmsgp;
#ifdef RESULTS_TRACKER
  Track *trkp;
#endif
#ifdef hashhit_dump_sortarray
  FILE *dumpfp;
  int dumpctr = 0;
#endif

  ERRMSG_CREATE(errmsgp);

  if (argc != 3) {
    printf("usage: %s <prefix hash index> <fastq reads> <FASTQ output file>\n",argv[0]);
    exit(0);
  }

  binfilnam = argv[1];
  readfilnam = argv[2];

#ifdef RESULTS_TRACKER
  if (NULL == (trkp = trackCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
#endif
#ifdef hashhit_dump_sortarray
  if (NULL == (dumpfp = EFOPEN(hashhit_helper_output_filnam, "wb")))
    ERRMSGNO(errmsgp, ERRCODE_NOFILE);
#endif

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  if (!(sqp = seqFastqCreate(0, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  printf("looking for read input file ...\n");
  sfp = seqIOopen(&errcode, readfilnam, SEQIO_READ, 0);
  if (errcode) ERRMSGNO(errmsgp, errcode);

  printf("Reading reference sequences ...\n");
  ssp = seqSetReadBinFil(&errcode, binfilnam);
  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  printf("Reading hash table ...\n");
  htp = hashTableRead(&errcode, binfilnam);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);
  hashTablePrintStats(stdout, htp);

  if (!(hhip = hashCreateHitInfo(HITINFO_BLKSZ, htp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(hhlp = hashCreateHitList(HASH_MAXN_HITS)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  sctr = 0;
  while(!seqIOstatus(sfp)) {
    printf("Processing read %i ...\n", sctr++);
    if ((errcode = seqFastqRead(sqp, sfp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = seqFastqEncode(sqp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    
    processRead(errmsgp,
		is_verbose,
#ifdef RESULTS_TRACKER
		trkp,
#endif
#ifdef hashhit_dump_sortarray
		dumpfp,
		&dumpctr,
#endif
		hhip, hhlp, htp, ssp, sqp, codecp);

#ifdef RESULTS_TRACKER
    if (!trackTestFlag(trkp, TRACKFLG_KMERHITS_OVERLAP)) {
      printf("%s\n", seqFastqGetSeqName(sqp));
    }
#endif
 }
  if (seqIOstatus(sfp) && 
      seqIOstatus(sfp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp));
  seqIOclose(sfp);

  hashDeleteHitList(hhlp);
  hashDeleteHitInfo(hhip);
  hashTableDelete(htp);
  seqFastqDelete(sqp);
  seqCodecDelete(codecp);
#ifdef RESULTS_TRACKER
  trackDelete(trkp);
#endif
#ifdef hashhit_dump_sortarray
  EFCLOSE(dumpfp);
#endif

  ERRMSG_END(errmsgp);
  exit(errcode);
}

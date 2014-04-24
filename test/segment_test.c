/** Testing of routines in segment.c
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

//#define segment_debug
//#define segment_test_verbose

#include "elib.h"
#include "sequence.h"
#include "hashidx.h"
#include "hashhit.h"
#include "segment.h"

enum {
  WITH_TERMCHAR = 1,
  REFSEQBLOCKSIZ = 64*1024*1024,
  HASH_MAXN_HITS = 128*1024*1024,
  HITINFO_BLKSZ =  1024,
  SEGLST_BLOCKSIZ = 10000,
  MASK32BIT = 0xFFFFFFFF,
  SHIFT_CUTOFF_NKTUP = 3,
  HASH_MIN_BASQUAL = 0,
};

static int checkSegmentList(SeqFastq *sqbufp, const SegLst *sglp, const SeqSet *ssp, 
			    const SeqFastq *sqp, const HashTable *htp, const SeqCodec *codecp)
{
  int errcode;
  unsigned char is_reverse;
  const char *qsp, *rsp;
  unsigned char nskip, ktup;
  int32_t len;
  uint32_t i, j, n_seeds, qo, qe, ss;
  SETSIZ_t so, se;
  
  ktup = hashTableGetKtupLen(htp, &nskip);
  is_reverse = segLstGetStats(sglp, NULL, &n_seeds, NULL);

  qsp = seqFastqGetConstSequence(sqp, NULL, NULL);
  for (i=0; i<n_seeds; i++) {
    len = segLstFetchSeed(&qo, &ss, i, sglp);
    if (is_reverse) {
      so = ss*nskip - len + ktup;
    } else {
      so *= nskip;
    }
#ifdef segment_test_verbose
    printf("[%i] qo = %u, so = %u, len = %i\n", i, qo, so, len);
#endif
    se = so + len;
    if ((errcode = seqSetFetchSegment(sqbufp, &so, &se, ssp, 
				      SEQCOD_ASCII)))
      return errcode;
    if (is_reverse &&
	(errcode = seqFastqReverse(sqbufp, codecp)))
      return errcode;
    qe = qo+len;
    rsp = seqFastqGetConstSequence(sqbufp, NULL, NULL);
    for (j = qo; j<qe && (*rsp); j++) {
      if (qsp[j] != *rsp++) {
	printf("Q:%s\n", qsp);
	printf("S:%s\n", seqFastqGetConstSequence(sqbufp, NULL, NULL));
	printf("%c:%c\n", qsp[j], *(rsp-1));
	break;
      }
    }
/*     printf("Q:%s\n", qsp); */
/*     printf("S:%s\n", seqFastqGetConstSequence(sqbufp, NULL, NULL)); */
    
    if (j<qe) {
      printf("SEED[%i] qo = %u, so = %llu, len = %i does not match at position %u ...", 
	     i, qo, (unsigned long long) so, len, j);
      return ERRCODE_FAILURE;
    }
  }
  return ERRCODE_SUCCESS;
}

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  unsigned char nskip;
  char *binfilnam, *readfilnam;
  int maxhit_per_tuple;
  unsigned int min_ktup;
  long sctr;
  SeqCodec *codecp;
  SeqFastq *sqp, *sqbufp;
  SeqIO *sfp;
  SeqSet *ssp;
  HashTable *htp;
  HashHitInfo *hhip;
  HashHitList *hhlp;
  SegLst *sglp;
#ifdef segment_with_segcand
  int max_depth, mincov_below_max;
  SegAliCands *sacp;
  unsigned char is_sensitive = 0;
#endif
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc != 7) {
    printf("usage: %s <hash table root name> <ncut> <nseeds>" \
	   " <mincov_below_max> <max_depth> <FASTQ file>\n",argv[0]);
    exit(0);
  }

  binfilnam = argv[1];
  maxhit_per_tuple = atoi(argv[2]);
  min_ktup = atoi(argv[3]);
#ifdef segment_with_segcand
   mincov_below_max = atoi(argv[4]);
  max_depth = atoi(argv[5]);
#endif
  readfilnam = argv[6];

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (!(sqbufp = seqFastqCreate(REFSEQBLOCKSIZ, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(sqp = seqFastqCreate(0, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(hhlp = hashCreateHitList(HASH_MAXN_HITS)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(sglp = segLstCreate(SEGLST_BLOCKSIZ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
#ifdef segment_with_segcand
  if (!(sacp = segAliCandsCreate(SEGLST_BLOCKSIZ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
#endif
  sfp = seqIOopen(&errcode, readfilnam, SEQIO_READ, 0);
  if (errcode) ERRMSGNO(errmsgp, errcode);

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

  if (!(hhip = hashCreateHitInfo(HITINFO_BLKSZ, htp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  printf("Working query sequence files ...\n");
  sctr = 0;
  while(!seqIOstatus(sfp)) {
#ifdef segment_with_segcand
    segAliCandsBlank(sacp);
#endif
    printf("Reading query sequence %li ...\n", sctr++);
    if ((errcode = seqFastqRead(sqp, sfp)))
      ERRMSGNO(errmsgp, errcode);
#ifdef segment_test_verbose
    printf("Collecting hits on forward strand ...\n");
#endif
    if ((errcode = seqFastqEncode(sqp, codecp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = hashCollectHitInfo(hhip, 0, HASH_MIN_BASQUAL, 0, 0, sqp, htp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = hashCollectHitsUsingCutoff(hhlp, maxhit_per_tuple, htp, hhip)))
      ERRMSGNO(errmsgp, errcode);

#ifdef segment_test_verbose
    hashPrintHitList(hhlp, stdout);
#endif
    if ((errcode = segLstFillHits(sglp, min_ktup, hhlp)))
      ERRMSGNO(errmsgp, errcode);
#ifdef segment_test_verbose
    segLstPrintSeeds(stdout, sglp);
#endif
    if ((errcode = seqFastqDecode(sqp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    if ((errcode = checkSegmentList(sqbufp, sglp, ssp, sqp, htp, codecp)))
      ERRMSGNO(errmsgp, errcode);
#ifdef segment_with_segcand
    if ((errcode = segAliCandsAddNoIndel(sacp, sglp)))
      ERRMSGNO(errmsgp, errcode);
#endif 
#ifdef segment_test_verbose 
    printf("Collecting hits on reverse complement strand ...\n");
#endif
    if ((errcode = seqFastqEncode(sqp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    
    if ((errcode = hashCollectHitInfo(hhip, 1, HASH_MIN_BASQUAL, 0, 0, sqp, htp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = hashCollectHitsUsingCutoff(hhlp, maxhit_per_tuple, htp, hhip)))
      ERRMSGNO(errmsgp, errcode);

#ifdef segment_test_verbose
    hashPrintHitList(hhlp, stdout);
#endif
    if ((errcode = segLstFillHits(sglp, min_ktup, hhlp)))
      ERRMSGNO(errmsgp, errcode);
#ifdef segment_test_verbose
    segLstPrintSeeds(stdout, sglp);
#endif
    if ((errcode = seqFastqDecode(sqp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    if ((errcode = checkSegmentList(sqbufp, sglp, ssp, sqp, htp, codecp)))
      ERRMSGNO(errmsgp, errcode);
#ifdef segment_with_segcand
    if ((errcode = segAliCandsAddNoIndel(sacp, sglp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = segAliCandsStats(sacp, mincov_below_max, max_depth, is_sensitive)))
      ERRMSGNO(errmsgp, errcode);
    segAliCandsPrint(stdout, max_depth, sacp);
#endif
  }
  if (seqIOstatus(sfp) && 
      seqIOstatus(sfp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp));
  seqIOclose(sfp);

  
  hashTableDelete(htp);
#ifdef segment_with_segcand
  segAliCandsDelete(sacp);
#endif
  segLstDelete(sglp);
  hashDeleteHitInfo(hhip);
  hashDeleteHitList(hhlp);
  seqSetDelete(ssp);
  seqFastqDelete(sqp);
  seqFastqDelete(sqbufp);
  seqCodecDelete(codecp);

  ERRMSG_END(errmsgp);

  return ERRCODE_SUCCESS;
}

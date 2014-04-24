/** Testing of alignment seeds lookup from hash table.
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
  HASH_MAXN_HITS = 1024*1024,
  HITINFO_BLKSZ =  1024,
  HASH_MIN_BASQUAL = 0,
};

typedef unsigned char UCHAR;

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  int sctr, maxhit_per_tuple;
  char is_verbose = 1, *reffilnam, *readfilnam;
  UCHAR seqset_flags = 0;
  char is_reverse = 0;
  UCHAR ktup, nskip, nbits_perf, nbits_key;
  SeqIO *sfp_in;
  SeqCodec *codecp;
  SeqFastq *sqp, *sqbufp;
  SeqSet *ssp;
  HashTable *htp;
  HashHitInfo *hhip;
  HashHitList *hhlp;
  ErrMsg *errmsgp;

  ERRMSG_CREATE(errmsgp);

  if (argc != 8) {
    printf("usage: %s <fasta subject> <ktup> <nskip> <nbits_key> <nbits_perf> <maxhit> <fastq reads>\n",argv[0]);
    exit(0);
  }
  
  reffilnam = argv[1];
  ktup = (UCHAR) atoi(argv[2]);
  nskip = (UCHAR) atoi(argv[3]);
  nbits_key = (UCHAR) atoi(argv[4]);
  nbits_perf = (UCHAR) atoi(argv[5]);
  maxhit_per_tuple = atoi(argv[6]);
  readfilnam = argv[7];

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (!(sqbufp = seqFastqCreate(SEQBLOCKSIZ, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  if (!(sqp = seqFastqCreate(0, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (!(ssp = seqSetCreate(SEQBLOCKSIZ, seqset_flags)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(htp = hashTableCreate(ktup, nskip, nbits_key, nbits_perf, 
			      (ktup*2 == nbits_key)? 
			      HASHIDXTYP_PERFECT: HASHIDXTYP_HASH32MIX)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(hhlp = hashCreateHitList(HASH_MAXN_HITS)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  sfp_in = seqIOopen(&errcode, reffilnam, SEQIO_READ, HASH_SEQREAD_BUFSIZ);
  if (errcode) ERRMSGNO(errmsgp, errcode);
  
  sctr = 0;
  while(!seqIOstatus(sfp_in)) {
    printf("Reading sequence %i ...\n", sctr++);
    if ((errcode = seqFastqRead(sqbufp, sfp_in)))
      ERRMSGNO(errmsgp, errcode);
    if ((errcode = seqFastqEncode(sqbufp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    if ((errcode = seqSetAddSequence(ssp, sqbufp)))
      ERRMSGNO(errmsgp, errcode);
  }
  if (seqIOstatus(sfp_in) && 
      seqIOstatus(sfp_in) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp_in));
  seqIOclose(sfp_in);

  printf("Setting up table ...\n");
  if ((errcode = hashTableSetUp(htp, sqbufp, ssp, NULL, codecp, NULL,
				is_verbose))) 
    ERRMSGNO(errmsgp, errcode);

  hashTablePrintStats(stdout, htp);

  if (!(hhip = hashCreateHitInfo(HITINFO_BLKSZ, htp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  /* work the reads */
  sfp_in = seqIOopen(&errcode, readfilnam, SEQIO_READ, 0);
  if (errcode) ERRMSGNO(errmsgp, errcode);
  sctr = 0;
  while(!seqIOstatus(sfp_in)) {
    printf("Reading read %i ...\n", sctr++);
    if ((errcode = seqFastqRead(sqp, sfp_in)))
      ERRMSGNO(errmsgp, errcode);
    if ((errcode = seqFastqEncode(sqp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    
    if ((errcode = hashCollectHitInfo(hhip, is_reverse, HASH_MIN_BASQUAL, 0, 0, sqp, htp)))
      ERRMSGNO(errmsgp, errcode);

    printf("Looking up seeds ...\n");
    hashBlankHitList(hhlp);
    if ((errcode = hashCollectHitsUsingCutoff(hhlp, maxhit_per_tuple, htp, hhip)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = hashCheckHitList(hhlp, sqp, htp, ssp, codecp)))
      ERRMSGNO(errmsgp, errcode);
  }
  if (seqIOstatus(sfp_in) && 
      seqIOstatus(sfp_in) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp_in));
  seqIOclose(sfp_in);

  hashDeleteHitList(hhlp);
  hashDeleteHitInfo(hhip);
  hashTableDelete(htp);
  seqSetDelete(ssp);
  seqFastqDelete(sqp);
  seqFastqDelete(sqbufp);
  seqCodecDelete(codecp);

  ERRMSG_END(errmsgp);
  exit(errcode);
}

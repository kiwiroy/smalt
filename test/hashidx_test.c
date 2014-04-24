/** Checking hash calculations
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
#include <string.h>

#include "elib.h"
#include "sequence.h"
#include "hashidx.h"

enum {
  HASH_SEQREAD_BUFSIZ = 128*1024*1024,
  SEQBLOCKSIZ = 64*1024*1024,
  NBITS_PERF = 8,
};

typedef unsigned char UCHAR;

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  char *hashfilnam, is_verbose = 1;
  UCHAR wordlen, nskip, nbits_key, hashtyp;
  UCHAR seqset_flags = SEQSET_COMPRESSED;
  SeqCodec *codecp;
  SeqFastq *sqbufp;
  SeqSet *ssp;
  HashTable *htp, *htBp;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc != 6) {
    printf("usage: %s <fasta subject> <wordlen> <nskip> <nbits_key> <binary file root name>\n", argv[0]);
    exit(0);
  }
  wordlen = (UCHAR) atoi(argv[2]);
  nskip = (UCHAR) atoi(argv[3]);
  nbits_key = (UCHAR) atoi(argv[4]);
  hashfilnam = argv[5];

  if (nbits_key >= wordlen*2) {
    printf("Construct a perfect hash index ...\n");
    hashtyp = HASHIDXTYP_PERFECT;
  } else {
    printf("Construct a hash index with collisions ...\n");
    hashtyp = HASHIDXTYP_HASH32MIX;
  }

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (!(sqbufp = seqFastqCreate(SEQBLOCKSIZ, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  if (!(ssp = seqSetCreate(SEQBLOCKSIZ, seqset_flags)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(htp = hashTableCreate(wordlen, nskip, nbits_key, NBITS_PERF, hashtyp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  printf("Reading reference sequences ...\n");
  if ((errcode = seqSetAddFromFastqFile(errmsgp, ssp, sqbufp, codecp, argv[1], is_verbose)))
    return errcode;

/*   printf("Compressing sequence set ...\n"); */
/*   if ((errcode = seqSetCompress(ssp, codecp))) */
/*     ERRMSGNO(errmsgp, errcode); */

  printf("Writing sequence set ...\n");
  if ((errcode = seqSetWriteBinFil(ssp, hashfilnam)))
    ERRMSGNO(errmsgp, errcode);

  printf("Setting up table ...\n");
  if ((errcode = hashTableSetUp(htp, sqbufp, ssp, NULL, codecp, 
				NULL, is_verbose))) 
    ERRMSGNO(errmsgp, errcode);

  hashTablePrintStats(stdout, htp);
  
  printf("Writing table to file ...\n");
  if ((errcode = hashTableWrite(hashfilnam, htp)))
    ERRMSGNO(errmsgp, errcode);

  printf("\nChecking table ...\n");
  errcode = hashTableCheckQuick(sqbufp, htp, ssp, codecp);
  if (errcode) {
    printf("Hash table broken.\n");
    ERRMSGNO(errmsgp, errcode);
  } else {
    printf("Hash table ok.\n");
  }
  
  printf("\nChecking table extensively ...\n");
  errcode = hashTableCheckExtensive(sqbufp, htp, ssp, codecp);
  if (errcode) {
    printf("Hash table broken.\n");
    ERRMSGNO(errmsgp, errcode);
  } else {
    printf("Hash table ok.\n");
  }
  
  printf("Reading table from file ... \n");
  htBp = hashTableRead(&errcode, hashfilnam);
  if (!(htBp) || (errcode))
    ERRMSGNO(errmsgp, errcode);

  printf("Comparing hash table ...\n");
  errcode = hashTableCmp(htp, htBp);
  if ((errcode)) {
    printf("Test failed.\n");
  } else {
    printf("Test ok.\n");
  }

  hashTableDelete(htBp);
  hashTableDelete(htp);
  seqSetDelete(ssp);
  seqFastqDelete(sqbufp);
  seqCodecDelete(codecp);

  ERRMSG_END(errmsgp);
  exit(errcode);
}

/** Testing of routines in segment.c
 *
 * Using segment.c::segAliCandsAddSimple() to add const-shift segments
 * to alignment candidates.
 */

#include <stdio.h>
#include <stdlib.h>

//#define segment_debug
#define segment_test_verbose

#include "elib.h"
#include "sequence.h"
#include "hashidx.h"
#include "hashhit.h"
#include "segment.h"

enum {
  REFSEQBLOCKSIZ = 64*1024*1024,
  HASH_MAXN_HITS = 128*1024*1024,
  HITINFO_BLKSZ = 1024,
  SEGLST_BLOCKSIZ = 10000,
  MAX_DEPTH = 500,
  HASH_MIN_BASQUAL = 0,
};

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  unsigned char ktup, nskip;
  char *binfilnam, *readfilnam;
  int maxhit_per_tuple;
  unsigned int min_ktup, min_cover;
  long sctr;
  SeqCodec *codecp;
  SeqFastq *sqp, *sqbufp;
  SeqIO *sfp;
  SeqSet *ssp;
  HashTable *htp;
  HashHitInfo *hhip;
  HashHitList *hhlp;
  SegLst *rev_sglp, *sglp;
  SegAliCands *sacp;
  SegQMask *sqmaskp;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc != 5) {
    printf("usage: %s <hash table root name> <ncut> <mincover>" \
	   " <FASTQ file>\n",argv[0]);
    exit(0);
  }

  binfilnam = argv[1];
  maxhit_per_tuple = atoi(argv[2]);
  min_cover = atoi(argv[3]);
  readfilnam = argv[4];

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

   if (!(rev_sglp = segLstCreate(SEGLST_BLOCKSIZ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  if (!(sqmaskp = segQMaskCreate(0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(sacp = segAliCandsCreate(SEGLST_BLOCKSIZ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

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
  ktup = hashTableGetKtupLen(htp, &nskip);

  printf("\nK-mer word length: %hi\n", (short) ktup);
  printf("Skip step size: %hi\n", (short) nskip);
  if (min_cover < ktup) min_cover = ktup;
  min_ktup = (min_cover-1)/ktup + 1;
  printf("Minimum cover: %u\n", min_cover);
  printf("Minimum number of k-tuples: %u\n\n", min_ktup);

  if (!(hhip = hashCreateHitInfo(HITINFO_BLKSZ, htp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  printf("Working query sequence files ...\n");
  sctr = 0;
  while(!seqIOstatus(sfp)) {
    segAliCandsBlank(sacp);

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

#ifdef segment_test_verbose 
    printf("Collecting hits on reverse complement strand ...\n");
#endif
    if ((errcode = hashCollectHitInfo(hhip, 1, HASH_MIN_BASQUAL, 0, 0, sqp, htp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = hashCollectHitsUsingCutoff(hhlp, maxhit_per_tuple, htp, hhip)))
      ERRMSGNO(errmsgp, errcode);

#ifdef segment_test_verbose
    hashPrintHitList(hhlp, stdout);
#endif
    if ((errcode = segLstFillHits(rev_sglp, min_ktup, hhlp)))
      ERRMSGNO(errmsgp, errcode);
#ifdef segment_test_verbose
    segLstPrintSeeds(stdout, rev_sglp);
#endif

    if ((errcode = segAliCandsAddFast(sacp, sqmaskp, sglp, min_cover, SEGCAND_UNKNOWN_SEQIDX)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = segAliCandsAddFast(sacp, sqmaskp, rev_sglp, min_cover, SEGCAND_UNKNOWN_SEQIDX)))
      ERRMSGNO(errmsgp, errcode);

#ifdef segment_test_verbose
    segAliCandsPrintRaw(stdout, MAX_DEPTH, sacp);
#endif
  }
  if (seqIOstatus(sfp) && 
      seqIOstatus(sfp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp));
  seqIOclose(sfp);

  
  hashTableDelete(htp);
  segAliCandsDelete(sacp);
  segQMaskDelete(sqmaskp);
  segLstDelete(rev_sglp);
  segLstDelete(sglp);
  hashDeleteHitList(hhlp);
  hashDeleteHitInfo(hhip);
  seqSetDelete(ssp);
  seqFastqDelete(sqp);
  seqFastqDelete(sqbufp);
  seqCodecDelete(codecp);

  ERRMSG_END(errmsgp);
 
  return ERRCODE_SUCCESS;
}

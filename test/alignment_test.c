/**< Test for full Smith-Waterman alignment routine 
 * in aligment.c 
 */

#include <stdlib.h>
#include <limits.h>

#include "../src/elib.h"
#include "../src/sequence.h"
#include "../src/score.h"
#include "../src/alignment.h"

enum {
  NLOOPS_TIMING = 10000,
};

int main(int argc, char *argv[])
{
  int errcode;
  int r_edge, l_edge;
  int q_start, q_end;
  int s_start, s_end;
  unsigned int uslen, qlen;
  char *infilnam;
  const char *psqp, *usqp;
  SeqIO *sfp;
  SeqFastq *seqAp, *seqBp;
  SeqCodec *codecp;
  ScorePenalties *spltyp;
  ScoreMatrix *smatp;
  ScoreProfile *sprofp;
  AliBuffer *bufp;
  AliRsltSet *rssp;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc < 8) {
    printf("usage: %s <fasta/fastq file [in]>" \
	   "<l_edge> <r_edge> <q_start> <q_end> <s_start> <s_end>\n",
	   argv[0]);
    exit(0);
  }

  infilnam = argv[1];
  l_edge  = atoi(argv[2]);
  r_edge  = atoi(argv[3]);
  q_start = atoi(argv[4]);
  q_end   = atoi(argv[5]);
  s_start = atoi(argv[6]);
  s_end   = atoi(argv[7]);

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(spltyp = scorePenaltiesCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(smatp = scoreCreateMatrix(codecp, spltyp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  
  if (!(sprofp = scoreCreateProfile(0, codecp, 
				    SCORPROF_SCALAR 
#ifdef SCORE_SIMD_SSE2
| SCORPROF_STRIPED_8 | SCORPROF_STRIPED_16
#endif
				    )))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(bufp = aliBufferCreate(0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(rssp = aliRsltSetCreate(NULL, 0, 0, 0, 0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(seqAp = seqFastqCreate(0, SEQTYP_FASTA)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(seqBp = seqFastqCreate(0, SEQTYP_FASTA)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);  

  sfp = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
  if (!sfp) ERRMSGNO(errmsgp, errcode);


  while(!(seqIOstatus(sfp))) {

    if ((errcode = seqFastqRead(seqAp, sfp)))
      ERRMSGNO(errmsgp, errcode);
    
    if ((errcode = seqFastqRead(seqBp, sfp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = seqFastqEncode(seqAp, codecp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = seqFastqEncode(seqBp, codecp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = scoreMakeProfileFromSequence(sprofp, seqBp, smatp)))
      ERRMSGNO(errmsgp, errcode);
    
    psqp = seqFastqGetConstSequence(seqBp, &qlen, NULL);
    usqp = seqFastqGetConstSequence(seqAp, &uslen, NULL);

    if ((errcode = aliBufferInit(bufp, qlen)))
      ERRMSGNO(errmsgp, errcode);

    if (uslen > INT_MAX)
      ERRMSGNO(errmsgp, ERRCODE_SEQLEN);

    if ((errcode = aliDebugFullSmiWat(rssp, bufp, sprofp, codecp, psqp, usqp, uslen,
				      -r_edge, -l_edge, 
				      s_start, s_end,
				      q_start, q_end)))
      ERRMSGNO(errmsgp, errcode);

  }

  if (seqIOstatus(sfp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp));

  seqIOclose(sfp);
  seqFastqDelete(seqBp);
  seqFastqDelete(seqAp);
  aliBufferDelete(bufp);
  aliRsltSetDelete(rssp);
  scoreDeleteProfile(sprofp);
  scoreDeleteMatrix(smatp);
  scorePenaltiesDelete(spltyp);
  seqCodecDelete(codecp);

  ERRMSG_END(errmsgp);

  return EXIT_SUCCESS;
}

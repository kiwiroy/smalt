/**< Test Smith-Waterman alignment routines 
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "elib.h"
#include "sequence.h"
#include "score.h"
#include "alibuffer.h"

int main(int argc, char *argv[])
{
  int errcode;
  int maxscor;
  long readctr;
  unsigned int uslen, qlen;
  char *fqfilA, *fqfilB;
  const char *psqp, *usqp;
  SeqIO *sfAp, *sfBp;
  SeqFastq *seqAp, *seqBp;
  SeqCodec *codecp;
  ScorePenalties *spltyp;
  ScoreMatrix *smatp;
  ScoreProfile *sprofp;
  AliBuffer *alibufp;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc != 3) {
    printf("usage: %s <sequence A [FASTA/FASTQ]> "\
	   "<sequence B [FASTA/FASTQ]>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  fqfilA = argv[1];
  fqfilB = argv[2];

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(spltyp = scorePenaltiesCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(smatp = scoreCreateMatrix(codecp, spltyp)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  
  if (!(sprofp = scoreCreateProfile(0, codecp, 
#ifdef SCORE_SIMD_SSE2
				    SCORPROF_STRIPED_8 | SCORPROF_STRIPED_16
#else
				    SCORPROF_STRIPED_32
#endif
				    )))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(alibufp = aliBufferCreate(0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(seqAp = seqFastqCreate(0, SEQTYP_FASTA)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(seqBp = seqFastqCreate(0, SEQTYP_FASTA)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);  

  if (NULL == (sfAp = seqIOopen(&errcode, fqfilA, SEQIO_READ, 0)))
    ERRMSGNO(errmsgp, errcode);
  
  if (NULL == (sfBp = seqIOopen(&errcode, fqfilB, SEQIO_READ, 0)))
    ERRMSGNO(errmsgp, errcode);
 
  while(!seqIOstatus(sfAp) && !seqIOstatus(sfBp)) {

    if ((errcode = seqFastqRead(seqAp, sfAp)))
      ERRMSGNO(errmsgp, errcode);
    
    if ((errcode = seqFastqRead(seqBp, sfBp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = seqFastqEncode(seqAp, codecp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = seqFastqEncode(seqBp, codecp)))
      ERRMSGNO(errmsgp, errcode);

    if ((errcode = scoreMakeProfileFromSequence(sprofp, seqBp, smatp)))
      ERRMSGNO(errmsgp, errcode);
    
    psqp = seqFastqGetConstSequence(seqBp, &qlen, NULL);
    usqp = seqFastqGetConstSequence(seqAp, &uslen, NULL);

    if ((errcode = aliBufferInit(alibufp, qlen)))
      ERRMSGNO(errmsgp, errcode);

    if (uslen > INT_MAX)
      ERRMSGNO(errmsgp, ERRCODE_SEQLEN);
 
    maxscor = 0;
    if ((errcode = swSIMDAlignStriped(&maxscor, alibufp, sprofp,
#if defined alignment_matrix_debug
				      codecp, psqp,
#endif	
				      usqp, uslen)))
      ERRMSGNO(errmsgp, errcode);

    printf("SW-score [%li]: %i\n", readctr, maxscor);
    readctr++;
  }

  if (seqIOstatus(sfAp) != ERRCODE_EOF && 
      seqIOstatus(sfAp) != ERRCODE_SUCCESS)
    ERRMSGNO(errmsgp, seqIOstatus(sfAp));
  if (seqIOstatus(sfBp) != ERRCODE_EOF &&
      seqIOstatus(sfBp) != ERRCODE_SUCCESS) 
    ERRMSGNO(errmsgp, seqIOstatus(sfBp));
 

  seqIOclose(sfBp);
  seqIOclose(sfAp);
  seqFastqDelete(seqBp);
  seqFastqDelete(seqAp);
  aliBufferDelete(alibufp);
  scoreDeleteProfile(sprofp);
  scoreDeleteMatrix(smatp);
  scorePenaltiesDelete(spltyp);
  seqCodecDelete(codecp);

  ERRMSG_END(errmsgp);

  return EXIT_SUCCESS;
}

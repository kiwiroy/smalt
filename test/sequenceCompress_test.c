/** Test sequence compression/decompression.
 * $Id: sequenceCompress_test.c,v 1.2 2009-05-29 16:45:49 hp3 Exp $
 *
 * Produce a reference file and a number of output files which should 
 * be identical to reference file (files are not compared)
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "elib.h"
#include "sequence.h"

#define MAXNAMLEN 256
#define LINEWIDTH 60

int main(int argc, char *argv[])
{
  char nambuf[MAXNAMLEN], qnambuf[MAXNAMLEN];
  int errcode = ERRCODE_SUCCESS;
  int s_start, s_end;
  long ctr = 0;
  SeqIO *sfp_in, *sfAp_out, *sfBp_out;
  SeqFastq *seqp, *seqBp;
  SeqCodec *codep;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 6) {
    printf("usage: %s <fasta/fastq file [in]> <fasta/fastq file 1 [out]> ",argv[0]);
    printf("<fasta/fastq file 2 [out]> <start> <end>\n");
    exit(0);
  }
  
  s_start = atoi(argv[4]);
  s_end = atoi(argv[5]);

  if (!(codep = seqCodecCreate()))
    ERRMSGNO(errmsg, ERRCODE_FAILURE);

  if (!(seqp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  if (!(seqBp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);
  
  sfp_in = seqIOopen(&errcode, argv[1], SEQIO_READ, 0);
  if (!sfp_in) ERRMSGNO(errmsg, errcode);
  
  sfAp_out = seqIOopen(&errcode, argv[2], SEQIO_WRITE_FASTQ, 0);
  if (!sfAp_out) ERRMSGNO(errmsg, errcode);

  sfBp_out = seqIOopen(&errcode, argv[3], SEQIO_WRITE_FASTQ, 0);
  if (!sfBp_out) ERRMSGNO(errmsg, errcode);  
   
  while(!(seqIOstatus(sfp_in) || 
	  seqIOstatus(sfAp_out) || 
	  seqIOstatus(sfBp_out))) {
    if ((errcode = seqFastqRead(seqp, sfp_in)))
      ERRMSGNO(errmsg, errcode);
    ctr++;
    //if ((ctr % 100) < 1) {
    //seqFastqGetSequence(seqp, &sl, NULL);
    //printf("sequence %d, length %d\n", ctr, sl);
    //}
    if ((errcode = seqFastqCheck(seqp)))
      ERRMSGNO(errmsg, errcode);

    strncpy(nambuf, seqFastqGetSeqName(seqp), MAXNAMLEN-1);
    nambuf[MAXNAMLEN-1] = '\0';
    if (seqFastqGetQualName(seqp)) {
      strncpy(qnambuf, seqFastqGetQualName(seqp), MAXNAMLEN-1);
      qnambuf[MAXNAMLEN-1] = '\0';
    } else {
      qnambuf[0] = '\0';
    }
    
    seqFastqBlank(seqBp);

    if ((errcode = seqFastqAppendSegment(seqBp, seqp, s_start, s_end, 0, NULL)))
      ERRMSGNO(errmsg, errcode);

    if ((errcode = seqFastqWrite(sfAp_out, seqBp, LINEWIDTH)))
      ERRMSGNO(errmsg, errcode);

    /* now encode/compress and decode */
    if ((errcode = seqFastqEncode(seqp, codep)))
      ERRMSGNO(errmsg, errcode);

    if ((errcode = seqFastqCheck(seqp)))
      ERRMSGNO(errmsg, errcode);

    if ((errcode = seqFastqCompress(seqp)))
      ERRMSGNO(errmsg, errcode);
   
    if ((errcode = seqFastqUncompress(seqBp, seqp, s_start, s_end, codep, 0)))
      ERRMSGNO(errmsg, errcode);

    if ((errcode = seqFastqWrite(sfBp_out, seqBp, LINEWIDTH)))
      ERRMSGNO(errmsg, errcode);
    
  }
  if (seqIOstatus(sfp_in) != ERRCODE_EOF) ERRMSGNO(errmsg, seqIOstatus(sfp_in));
  if (seqIOstatus(sfAp_out)) ERRMSGNO(errmsg, seqIOstatus(sfAp_out));
  if (seqIOstatus(sfBp_out)) ERRMSGNO(errmsg, seqIOstatus(sfBp_out));

  seqIOclose(sfBp_out);
  seqIOclose(sfAp_out);
  seqIOclose(sfp_in);
  seqFastqDelete(seqBp);
  seqFastqDelete(seqp);
  seqCodecDelete(codep);

  ERRMSG_END(errmsg);

  return ERRCODE_SUCCESS;
}
  

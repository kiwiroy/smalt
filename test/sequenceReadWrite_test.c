/** Test routine for Read/Write of Fastq files
 */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Genome Research Ltd.                                  * 
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
#include <string.h>
#include <stdlib.h>
#include "elib.h"
#include "sequence.h"

#define MAXNAMLEN 256
#define OUTPUT_LINEWIDTH 60

int main(int argc, char *argv[])
{
  char nambuf[MAXNAMLEN], qnambuf[MAXNAMLEN];
  const char *cp;
  int errcode = ERRCODE_SUCCESS;
  int encode_flag;
  long ctr = 0;
  unsigned int sl;
  SeqIO *sfp_in, *sfp_out;
  SeqFastq *seqp, *seqBp;
  SeqCodec *codep;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 4) {
    printf("usage: %s <fasta/fastq file [in]> <fasta/fastq file [out]> <encode flag>\n",argv[0]);
    exit(0);
  }
  
  encode_flag = atoi(argv[3]);

  if (!(codep = seqCodecCreate()))
    ERRMSGNO(errmsg, ERRCODE_FAILURE);

  if (!(seqp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  if (!(seqBp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);
  

  sfp_in = seqIOopen(&errcode, argv[1], SEQIO_READ, 0);
  if (!sfp_in) ERRMSGNO(errmsg, errcode);
  
  sfp_out = seqIOopen(&errcode, argv[2], SEQIO_WRITE_FASTQ, 0);
  if (!sfp_out) ERRMSGNO(errmsg, errcode);
  
  while(!(seqIOstatus(sfp_in) || seqIOstatus(sfp_out))) {
    if ((errcode = seqFastqRead(seqp, sfp_in)))
      ERRMSGNO(errmsg, errcode);
    ctr++;
    //if ((ctr % 100) < 1) {
    seqFastqGetSequence(seqp, &sl, NULL);
    //printf("sequence %d, length %d\n", ctr, sl);
    //}
    if ((errcode = seqFastqCheck(seqp)))
      ERRMSGNO(errmsg, errcode);
    if (encode_flag) {
      strncpy(nambuf, seqFastqGetSeqName(seqp), MAXNAMLEN-1);
      nambuf[MAXNAMLEN-1] = '\0';
      if ((cp = seqFastqGetQualName(seqp))) {
	strncpy(qnambuf, seqFastqGetQualName(seqp), MAXNAMLEN-1);
	qnambuf[MAXNAMLEN-1] = '\0';
      } else {
	qnambuf[0] = '\0';
      }
      if ((errcode = seqFastqEncode(seqp, codep)))
	ERRMSGNO(errmsg, errcode);

      if ((errcode = seqFastqCheck(seqp)))
	ERRMSGNO(errmsg, errcode);

      if (encode_flag == 2) {
	if ((errcode = seqFastqCompress(seqp)))
	  ERRMSGNO(errmsg, errcode);
	
	if ((errcode = seqFastqUncompress(seqBp, seqp, 0, 0, codep, 0)))
	  ERRMSGNO(errmsg, errcode);
	
	seqFastqBlank(seqp);
	
	if ((errcode = seqFastqAppendSegment(seqp, seqBp, 0, 0, 0, codep)))
	  ERRMSGNO(errmsg, errcode);

	if ((errcode = seqFastqSetAscii(seqp, nambuf, NULL, qnambuf, NULL)))
	  ERRMSGNO(errmsg, errcode);
      } else {
	if ((errcode = seqFastqDecode(seqp, codep)))
	  ERRMSGNO(errmsg, errcode);    
      }
      if ((errcode = seqFastqCheck(seqp)))
	ERRMSGNO(errmsg, errcode);  
    }
    if ((errcode = seqFastqWrite(sfp_out, seqp, OUTPUT_LINEWIDTH)))
      ERRMSGNO(errmsg, errcode);
  }
  if (seqIOstatus(sfp_in) != ERRCODE_EOF) ERRMSGNO(errmsg, seqIOstatus(sfp_in));
  if (seqIOstatus(sfp_out)) ERRMSGNO(errmsg, seqIOstatus(sfp_out));

  seqIOclose(sfp_out);
  seqIOclose(sfp_in);
  seqFastqDelete(seqBp);
  seqFastqDelete(seqp);
  seqCodecDelete(codep);

  ERRMSG_END(errmsg);

  return ERRCODE_SUCCESS;
}
  

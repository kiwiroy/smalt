/** Test German's bambam library for reading read pairs from
 * SAM and bam files.
 */

#include <bambamc/BamBam_BamCollatorInterface.h>
#include <stdio.h>
#include <stdlib.h>
#include "elib.h"

static char const filtyp[] = "bam"; /* can be "bam" or "sam" */
static char const tmpdirnam[] = "tmpdir_bamcoll"; /* directory for temporary ouput files */
static int const keep_orphans = 1;

int main (int argc, char *argv[])
{
  int id;
  ErrMsg *errmsgp;
  BamBam_FastQRead readA, readB;

  ERRMSG_CREATE(errmsgp);

  if (argc < 2) {
    printf("usage: %s <bam file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  char *bamfilnam = argv[1];


  id = BamBam_AllocBamCollator(bamfilnam, filtyp, tmpdirnam, keep_orphans);
  
  if (id < 0) 
    ERRMSGNO(errmsgp, ERRCODE_NOFILE);

  while ( BAMBAM_ALIGNMENT_TYPE_NONE != BamBam_ReadPair(id, &readA, &readB, 
							NULL, NULL, '\n') ) {
    printf("New pair\n");
    printf("readA: %s\n", (readA.name)? readA.name: "None");
    printf("readB: %s\n", (readB.name)? readB.name: "None"); 
  }
  
  BamBam_FreeBamCollator(id);

  ERRMSG_END(errmsgp);

  exit(EXIT_SUCCESS);
}

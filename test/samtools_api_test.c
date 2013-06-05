/** Test SAMTOOLS API
 */

#include <sam.h>
#include <stdio.h>
#include <stdlib.h>
#include "elib.h"

int main (int argc, char *argv[])
{
  samfile_t *samfp;
  char const *bamfilnam;
  bam1_t *alip;
  ErrMsg *errmsgp;

  ERRMSG_CREATE(errmsgp);

  if (argc < 2) {
    printf("usage: %s <bam file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  bamfilnam = argv[1];

  EMALLOCP0(alip);
  if (alip == NULL)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  samfp = samopen(bamfilnam, "rb", NULL);
  
  
  samclose(samfp);

  ERRMSG_END(errmsgp);
  exit(EXIT_SUCCESS);
}

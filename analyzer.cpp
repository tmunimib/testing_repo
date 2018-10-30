#include "kseq.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

/*****************************************
 * Main
 *****************************************/

int main(int argc, char *argv[]) {
  gzFile id_file;
  kseq_t *seq;
int file_line;

id_file =
      gzopen("example/", "r");
  seq = kseq_init(transcript_file);

  // open and read the .fa
  while ((file_line = kseq_read(seq)) >= 0) { 


	}
	printf("return value: %d\n", file_line);
  kseq_destroy(seq);
  gzclose(id_file);
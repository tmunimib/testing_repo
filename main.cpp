#include "bloomfilter.h"
#include "kseq.h"
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector.hpp" // for the bit_vector class
#include "sdsl/util.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <seqan/align.h>
#include <stdio.h>
#include <string>
#include <vector> // std::vector
#include <zlib.h>

const size_t sizebloom = 1000000;
BF bloom(sizebloom);

// function search, returns indexes in a vector
vector<int> search(const string &kmer, BF &bloomfilter) {
  return bloomfilter.get_index(kmer);
}
using namespace seqan;

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

/*****************************************
 * Main
 *****************************************/

int main(int argc, char *argv[]) {
  gzFile transcript_file;
  gzFile read_file;
  kseq_t *seq_tr;
  kseq_t *seq_re;
  int file_line1;
  int file_line2;
  map<int, string> legend_ID;
  int mapped_ID = 0;
  const int kmer_length = 60;
  vector<string> transcript_kmers_vec;
  BF bloom(sizebloom);
  string name_transcript;
  string transcript_name = "";
  string read_name = "";
  typedef String<char> TSequence;             // sequence type
  typedef Align<TSequence, ArrayGaps> TAlign; // align type

  if (argc > 1) {
    transcript_name = argv[1];
    read_name = argv[2];
  } else {
    cout << "Error in input" << endl;
    return 1;
  }

  transcript_file =
      gzopen(transcript_name.c_str(), "r"); // STEP 2: open the file handler
  seq_tr = kseq_init(transcript_file);      // STEP 3: initialize seq

  // open and read the .fa
  while ((file_line1 = kseq_read(seq_tr)) >=
         0) { // STEP 4: read sequence of transcript
    string transcript_string;

    read_file = gzopen(read_name.c_str(), "r"); // STEP 2: open the file handler
    seq_re = kseq_init(read_file);              // STEP 3: initialize seq

    // open and read the .fa
    while ((file_line2 = kseq_read(seq_re)) >=
           0) { // STEP 4: read sequence of transcript
      string read_string;

      ifstream file;
      file.open("mismatch.txt");
      std::string line;
      std::string partial;

      std::vector<std::string> tokens;

      while (std::getline(file, line)) { // '\n' is the default delimiter

        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token,
                            '\t')) { // but we can specify a different one
          tokens.push_back(token);
          cout << token << endl;
        }
      }
    }
  }

  printf("return value: %d\n", file_line1);
  kseq_destroy(seq_tr);     // STEP 5: destroy seq
  gzclose(transcript_file); // STEP 6: close the file handler
  printf("return value: %d\n", file_line2);
  kseq_destroy(seq_re); // STEP 5: destroy seq
  gzclose(read_file);   // STEP 6: close the file handler

  /****************************************
   *ALIGNMENT
   ****************************************/

  /*
    // read fasta created as output in the preious step
    gzFile test_file;

    test_file =
  <<<<<<< HEAD
        gzopen("final_id.fa", "r"); // STEP 2: open the file handler
  =======
        gzopen("example/final_id.fa", "r"); // STEP 2: open the file handler
  >>>>>>> 5576e7ea74b2493708191191bea3ca6c0bd994c4
    seq = kseq_init(test_file);             // STEP 3: initialize seq

    // seq1--> read
    // seq2--> transcript

    ofstream al_file;
    al_file.open("alignment.fa");
    while ((file_line = kseq_read(seq)) >= 0) { // read final_id
      kseq_t *seq_tran;
      string transcript_assigned;
      transcript_assigned = seq->name.s;
      TSequence seq1 = seq->seq.s;
      gzFile transcript_global;
      string name_tr;
          TSequence seq2;

      transcript_global = gzopen("example/chrY_mod.fa", "r");
      seq_tran = kseq_init(transcript_global);
      while ((file_line = kseq_read(seq_tran)) >= 0) { // read final_id
        name_tr = seq_tran->name.s;
        if (name_tr == transcript_assigned) {
          seq2 = seq_tran->seq.s;
          Align<String<char>> ali;
          resize(rows(ali), 2);
          assignSource(row(ali, 0), seq1);
          assignSource(row(ali, 1), seq2);
          al_file << ">" << name_tr << endl;
          al_file << ";Score = "
                  << localAlignment(ali, Score<int>(3, -3, -2, -2),
  DynamicGaps())
  <<<<<<< HEAD
                  << ";" << endl;
  =======
                  << endl;
  >>>>>>> 5576e7ea74b2493708191191bea3ca6c0bd994c4
          al_file << seq1 << endl;
        }
      }
      kseq_destroy(seq_tran);     // STEP 5: destroy seq
      gzclose(transcript_global); // STEP 6: close the file handler
    }
    al_file.close();
    printf("return value: %d\n", file_line);
    kseq_destroy(seq);  // STEP 5: destroy seq
    gzclose(test_file); // STEP 6: close the file handler
  <<<<<<< HEAD


    //read fasta with reads and scores of their alignment with transcript
    gzFile align_fasta;

    align_fasta =
        gzopen("example/alignment.fa", "r");
    seq = kseq_init(align_fasta);
    while ((file_line = kseq_read(seq)) >= 0) {
  =======
  */
  /*
    // read fasta with reads and scores of their alignment with transcripts
    gzFile align_fasta;
  kseq_t *seq;
    align_fasta = gzopen("alignment.fa", "r"); // STEP 2: open the file handler
    seq = kseq_init(align_fasta);
    while ((file_line1 = kseq_read(seq)) >= 0) { // read final_id

      cout << seq->name.s << endl;
      cout << seq->comment.l << endl;
      cout << "seq: " << seq->seq.s << endl;
    }
    printf("return value: %d\n", file_line1);

    kseq_destroy(seq);
    gzclose(align_fasta);


        kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(align_fasta);  // STEP 6: close the file handler
  */
  return 0;
}

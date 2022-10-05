//-----------------------------------------------------------------------------
// "Manually" calculate metric values directly from original SNP (tped) file.
// This uses The CoMet output metrics file ONLY to get the required metrics
// coords, not the actual metrics values.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <string>
#include <sstream>
#include <iostream>

enum {NUM_WAY_MAX = 3};

#define ASSERT(condition, command) { if (!(condition)) {command; return 1;} }

enum{MAX_LABEL_LEN = 255};

//-----------------------------------------------------------------------------

bool is_delim(int c) {
  return '\t' == c || ' ' == c;
}

//-----------------------------------------------------------------------------

template<class T>
T min(T x, T y) {
  return x < y ? x : y;
}

//-----------------------------------------------------------------------------

template<class Field_t, class Sum_t>
int calculate_metric_elt(char* metric_type_prec, int num_way,
  int vector_num[NUM_WAY_MAX], int bit_num[NUM_WAY_MAX],
  FILE* snptxtfile, FILE* lifile) {

  // Initializations.

  const bool is_ccc = strcmp(metric_type_prec, "ccc") == 0;
  const bool is_duo = strcmp(metric_type_prec, "duo") == 0;
  const bool is_czek = strcmp(metric_type_prec, "czekanowski_single") == 0 ||
                       strcmp(metric_type_prec, "czekanowski_double") == 0;
  const bool is_single = strcmp(metric_type_prec, "czekanowski_single") == 0;

  //CHECK
  // Allele labels would be "A", "C", "T", "G", and if none then use "X".
  // NOTE: we will use the term "allele label" to denote a single one of these letters.

  const int allele_label_null = 'X';

  // Initialize line labels to null.



// TODO: use ASSERT everywhere.


// TODO: replace with c++ string.



  const int line_label_len = MAX_LABEL_LEN;
  unsigned char line_label[NUM_WAY_MAX][line_label_len+1];
  for (int i=0; i<num_way; ++i)
    for (int j=0; j<line_label_len+1; ++j)
      line_label[i][j] = 0;

  // Access the snp txt file.

  int fseek_success = fseek(snptxtfile, 0, SEEK_SET);
  ASSERT(0 == fseek_success,
         fprintf(stderr, "Error: error reading SNP data file (0).\n"));

  // Allele label is represented by a character, stored as an int.

  Field_t c = 0;

  const int num_allele_labels_per_field = is_ccc ? 2 : 1;

  // Allocate storage for allele labels and allele bits for each needed SNP.

  enum {NUM_ALLELE_LABELS_PER_FIELD_MAX = 2};

  Field_t allele_labels[NUM_WAY_MAX][NUM_ALLELE_LABELS_PER_FIELD_MAX];
  Field_t* elts[NUM_WAY_MAX];
  for (int way_num=0; way_num<num_way; ++way_num) {
    allele_labels[way_num][0] = allele_label_null;
    allele_labels[way_num][1] = allele_label_null;
    elts[way_num] = NULL;
  }

  int num_field = 0;

  // Initial size of snp file line - will be reallocated as needed.

  size_t linesize = 1;
  unsigned char* line = (unsigned char*)malloc(linesize);

  //----------
  // LOOP over num_way to input the required vectors.

  for (int way_num=0; way_num<num_way; ++way_num) {

    // Get location in the line index file needed to get the line index.

    fseek_success = fseek(lifile, vector_num[way_num]*sizeof(size_t), SEEK_SET);
    ASSERT(0 == fseek_success,
           fprintf(stderr, "Error: error reading line_indices file (1).\n"));

    // Get the index to the required line in the snptxtfile.

    size_t line_index = 0;
    size_t num_read = fread(&line_index, sizeof(size_t), 1, lifile);
    ASSERT(1 == num_read,
           fprintf(stderr, "Error: error reading line_indices file (2).\n"));

    // Set next-read pointer to start of needed line in snptxtfile.

    fseek_success = fseek(snptxtfile, line_index, SEEK_SET);
    ASSERT(0 == fseek_success,
           fprintf(stderr, "Error: error reading SNP data file (1).\n"));

    const int num_frontmatter_fields = 4;

    // Read the line from snptxtfile - read up to newline.

    const size_t line_len = getline((char**)&line, &linesize, snptxtfile) - 1;
    ASSERT(line && line_len <= linesize,
      fprintf(stderr, "Error: memory %zu %zu\n", line_len, linesize));

    // Allocate elts[way_num] since now we have an upper bound on size.
    elts[way_num] = (Field_t*)malloc(line_len *
      num_allele_labels_per_field * sizeof(Field_t));

    //----------
    // LOOP over chars in line.
    //----------

    for (size_t i=0, char_index=0, num_delim=0; i<line_len; ++i) {

      c = line[i];

      // Skip tab or space (these are separators).

      if (is_delim(c)) {
        ++num_delim;
        continue;
      }

      // Get line label if at column 2.

      if (num_delim == 1) {
        // Append character to label.
        line_label[way_num][char_index++] = c; //check
        continue;
      }

      // If finished with frontmatter then reset index.

      if (num_delim == num_frontmatter_fields - 1)
        char_index = 0;

      // Finished processing frontmatter.

      if (num_delim < num_frontmatter_fields)
        continue;

      //----------
      // 1. Handle Czekanowski case.
      //----------

      if (is_czek) {
        // Find end of token.
        int len = 1;
        while (i+len < line_len && !is_delim(line[i+len])) {
          ++len;
        }

        std::string line_((char*)line);

        // Write token to file in binary format.
        std::string token_string(line_.substr(i, len));
        std::stringstream token_stream(token_string);
        Field_t token;
        token_stream >> elts[way_num][char_index++];
        if (0 == way_num)
          ++num_field;

//printf("%i %i %f\n", way_num, char_index-1, (double) elts[way_num][char_index-1]);

        // Complete the processing of this token.
        i += len - 1;

//printf("%i %f   %i\n", (int)vector_num[way_num], (float)elts[way_num][char_index-1], (int)(char_index-1));


      } else {

      //----------
      // 2. Handle ccc, duo case.
      //----------

        // Store allele bit.

        elts[way_num][char_index++] = c;

        // If first cycle of num_way loop, then count num_field.

        if (0 == way_num && char_index % num_allele_labels_per_field == 0)
          ++num_field;

        // Record allele label. Normalize into alphabetical order.

        if (c != '0') {
          if (allele_labels[way_num][0] == allele_label_null) {
            // Store first of 2.
            allele_labels[way_num][0] = c;
          } else if (allele_labels[way_num][1] == allele_label_null) {
            // Store second of 2.
            if (allele_labels[way_num][0] < c) {
              // No alpha sort.
              allele_labels[way_num][1] = c;
            } else if (allele_labels[way_num][0] > c) {
              // Alpha sort.
              allele_labels[way_num][1] = allele_labels[way_num][0];
              allele_labels[way_num][0] = c;
            }
          }
        } // if c

      } // if (is_czek)

    } // for i

  } // for way_num

  free(line);

  // We now have all 2 (or 3) vectors and their associated data.

  // Now compute metric value.

  //----------
  // LOOP over way_num to get sum_i's and count_i's
  //----------

  Sum_t count1[NUM_WAY_MAX] = {};
  Sum_t sum1[NUM_WAY_MAX] = {};

  for (int way_num=0; way_num<num_way; ++way_num) {
    count1[way_num] = 0;
    sum1[way_num] = 0;

    //----------
    if (is_czek) {
    //----------

      count1[way_num] = num_field;
      sum1[way_num] = 0;

      #pragma omp parallel for reduction(+:sum1[way_num])
      for (int f=0; f<num_field; ++f) {
        sum1[way_num] += elts[way_num][f];
//printf("  %i %i %f %f\n", way_num, f, (float)elts[way_num][f], (float)sum1[way_num]);
      }

    //----------
    } else if (is_duo) {
    //----------

      #pragma omp parallel for reduction(+:sum1[way_num]) reduction(+:count1[way_num])
      for (int f=0; f<num_field; ++f) {
        const Field_t e0 = elts[way_num][f];

        if (e0 == '0')
          continue;

        count1[way_num] += 1;

        const int rho = (e0 == allele_labels[way_num][bit_num[way_num]]);
        sum1[way_num] += rho;
      } // for f

    //----------
    } else { // if (!is_duo && !is_czek) // if (is_ccc)
    //----------

      #pragma omp parallel for reduction(+:sum1[way_num]) reduction(+:count1[way_num])
      for (int f=0; f<num_field; ++f) {
        const Field_t e0 = elts[way_num][2*f];
        const Field_t e1 = elts[way_num][2*f+1];

        if (e0 == '0')
          continue;

        count1[way_num] += 1;

        // Calculate row and add to sum - see paper

        const int rho = (e0 == allele_labels[way_num][bit_num[way_num]]) +
                        (e1 == allele_labels[way_num][bit_num[way_num]]);
        sum1[way_num] += rho;
      } // for f

    //----------
    } // if (is_czek)
    //----------

  } // for way_num

  // Now get sum_{ij}'s (or sum_{ijk}'s if 3-way)

  Sum_t countijk = 0;
  Sum_t sumijk = 0;

  //----------
  if (is_czek) {
  //----------

    #pragma omp parallel for reduction(+:sumijk)
    for (int f=0; f<num_field; ++f) {

      const Sum_t e0 = elts[0][f];
      const Sum_t e1 = elts[1][f];
      const Sum_t e2 = 3 == num_way ? elts[2][f] : (Sum_t)0;

      sumijk += 2 == num_way ? min(e0, e1) :
        min(e0, e1) + min(e1, e2) + min(e2, e0) - min(e0, min(e1, e2));;

    } // for f

  //----------
  } else if (is_duo) {
  //----------

    #pragma omp parallel for reduction(+:sumijk) reduction(+:countijk)
    for (int f=0; f<num_field; ++f) {
      const int e00 = elts[0][f];
      if (e00 == '0')
        continue;

      const int e10 = elts[1][f];
      if (e10 == '0')
        continue;

      const int e20 = num_way == 2 ? 1 : elts[2][f];
      if (num_way == 3 && e20 == '0')
        continue;

      countijk += 1;

      const int rho0 = (e00 == allele_labels[0][bit_num[0]]);
      const int rho1 = (e10 == allele_labels[1][bit_num[1]]);
      const int rho2 = num_way == 2 ? 1 :
                       (e20 == allele_labels[2][bit_num[2]]);

      sumijk += rho0 * rho1 * rho2;
    } // for f

  //----------
  } else { // if (!is_duo && !is_czek) // if (is_ccc)
  //----------

    #pragma omp parallel for reduction(+:sumijk) reduction(+:countijk)
    for (int f=0; f<num_field; ++f) {
      const int e00 = elts[0][2*f];
      const int e01 = elts[0][2*f+1];
      if (e00 == '0')
        continue;

      const int e10 = elts[1][2*f];
      const int e11 = elts[1][2*f+1];
      if (e10 == '0')
        continue;

      const int e20 = num_way == 2 ? 1 : elts[2][2*f];
      const int e21 = num_way == 2 ? 1 : elts[2][2*f+1];
      if (num_way == 3 && e20 == '0')
        continue;

      countijk += 1;

      const int rho0 = (e00 == allele_labels[0][bit_num[0]]) +
                       (e01 == allele_labels[0][bit_num[0]]);
      const int rho1 = (e10 == allele_labels[1][bit_num[1]]) +
                       (e11 == allele_labels[1][bit_num[1]]);
      const int rho2 = num_way == 2 ? 1 :
                       (e20 == allele_labels[2][bit_num[2]]) +
                       (e21 == allele_labels[2][bit_num[2]]);

      sumijk += rho0 * rho1 * rho2;
    } // for f

  //----------
  } // if (is_czek)
  //----------

  // substitute into formula

  double value = 0;

  //----------
  if (is_czek) {
  //----------

    value = (2 == num_way ? 2 : 3./2) * sumijk / (sum1[0] + sum1[1] + sum1[2]);

//printf("%i %i %i %f %f %f %f\n", (int)vector_num[0], (int)vector_num[1], (int)vector_num[2], (float)sumijk, (float)sum1[0], (float)sum1[1], (float)sum1[2]);

  } else { // if (!is_czek)

    const int cbpe = is_duo ? 1 : 2; // counted bits per element.
    const int cbpe_n = 2 == num_way ? cbpe * cbpe : cbpe * cbpe * cbpe;

    double f1[NUM_WAY_MAX];
    for (int way_num=0; way_num<num_way; ++way_num)
      f1[way_num] = 0 == count1[way_num] ? 0 :
        sum1[way_num] * 1. / (double)( cbpe * count1[way_num] );

    const double fijk = 0 == countijk ? 0 :
                     sumijk * 1. / (double)( cbpe_n * countijk );

    // NOTE hard-wired constant here

    const double multiplier = is_duo ? (double)4. : 9. / (double)2.;
    const double param = 2. / (double)3.;

    value = 2 == num_way ?
      multiplier * fijk * (1 - param * f1[0]) *  (1 - param * f1[1]) :
      multiplier * fijk * (1 - param * f1[0]) *  (1 - param * f1[1])
                        * (1 - param * f1[2]);

  //----------
  } // if (is_czek)
  //----------

  // Permute labels to output each result with a uniform order of labels.
  // By convention assume output line nums increasing, e.g. "0 1" not "1 0".

  struct Perm {
    Perm(int v0, int v1) : data{v0, v1, 0} {}
    Perm(int v0, int v1, int v2) : data{v0, v1, v2} {}
    int& operator[](int i) {return data[i];}
    int data[NUM_WAY_MAX];
  };

  Perm perm = 2 == num_way
    ? (vector_num[0] > vector_num[1] ?  Perm(0, 1) :
                                        Perm(1, 0))
    : (vector_num[0] > vector_num[1] && 
       vector_num[1] > vector_num[2] ?  Perm(0, 1, 2) :
       vector_num[0] > vector_num[2] && 
       vector_num[2] > vector_num[1] ?  Perm(0, 2, 1) :
       vector_num[1] > vector_num[0] && 
       vector_num[0] > vector_num[2] ?  Perm(1, 0, 2) :
       vector_num[1] > vector_num[2] && 
       vector_num[2] > vector_num[0] ?  Perm(1, 2, 0) :
       vector_num[2] > vector_num[0] && 
       vector_num[0] > vector_num[1] ?  Perm(2, 0, 1) :
                                        Perm(2, 1, 0));

  // Do output to stdout

  for (int way_num=0; way_num<num_way; ++way_num) {
    int iperm = perm[way_num];
    printf(0 != way_num ? " " : "");
    printf("%i", vector_num[iperm]);
    if (!is_czek)
      printf(" %i", bit_num[iperm]);
  } // for way_num

  for (int way_num=0; way_num<num_way; ++way_num) {
    int iperm = perm[way_num];
    printf(" %s", line_label[iperm]);
    if (!is_czek)
      printf("_%c", allele_labels[iperm][bit_num[iperm]]);
  } // for way_num

  printf(" %f\n", value);

  for (int way_num=0; way_num<num_way; ++way_num) {
    free(elts[way_num]);
  }
  return 0;
}

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  // Help message.

  if (argc < 4) {
    printf("validate: create metrics data for validation of calculations\n");
    printf("Usage: validate "
           "<metric_type_pre> <num_way> <snptxtfile> <line_index_file>\n");
    printf("Here stdin has the (ascii) metric entries, one per line.\n");
    return 0;
  }

  // Parse arguments.

  int argnum = 1;
  char* metric_type_prec = argv[argnum];

  if (strcmp(metric_type_prec, "ccc") != 0 &&
      strcmp(metric_type_prec, "duo") != 0 &&
      strcmp(metric_type_prec, "czekanowski_single") != 0 &&
      strcmp(metric_type_prec, "czekanowski_double") != 0) {
    fprintf(stderr, "Error: invalid metric_type_prec\n");
    return 1;
  } 

  const bool is_czek = strcmp(metric_type_prec, "czekanowski_single") == 0 ||
                       strcmp(metric_type_prec, "czekanowski_double") == 0;
  const bool is_single = strcmp(metric_type_prec, "czekanowski_single") == 0;
  ++argnum;

  if (strcmp(argv[argnum], "2") != 0 && strcmp(argv[argnum], "3") != 0) {
    fprintf(stderr, "Error: invalid num_way\n");
    return 1;
  }
  const int num_way = atoi(argv[argnum++]);

  FILE* snptxtfile = fopen(argv[argnum], "r");
  if (!snptxtfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  ++argnum;

  FILE* lifile = fopen(argv[argnum], "rb");
  if (!lifile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  ++argnum;

  // Prepare to read from stdin.

  enum {LINE_LEN_MAX = 4096};
  unsigned char line[LINE_LEN_MAX];
  //enum {LINE_LEN_MAX = 2000000};
  //unsigned char* line = (unsigned char*)malloc(LINE_LEN_MAX * sizeof(*line));
  //if (!line) {
  //  fprintf(stderr, "Error in malloc");
  //  return 1;
  //}

  int lineptr = 0;
  int argc_ = 1;
  char* argv_[LINE_LEN_MAX];
  int c = 0;

  // Loop over chars from stdin

  while ((c = fgetc(stdin)) != EOF) {

    // If not newline then append char to line.

    if (c != '\n') {
      line[lineptr++] = c;
      if (lineptr >= LINE_LEN_MAX) {
        fprintf(stderr, "Error: line too long");
        return 1;
      }
      continue;
    }

    // Have entire line, now create an (argc, argv) arg list consisting of the tokens.

    line[lineptr] = 0;
    argv_[0] = (char*)&line[lineptr];
    for (int i=0; i<lineptr; ++i) {
      if (is_delim(line[i]))
        line[i] = 0;
    }
    for (int i=0; i<lineptr; ++i) {
      if (line[i] != 0 && (i == 0 || line[i-1] == 0))
          argv_[argc_++] = (char*)&line[i];
    }
    argc_ = argc_ < 2*num_way+1 ? argc_ : 2*num_way+1;

    // Process this line.

    // Only use the first several tokens, specifying line and bit numbers;
    // discard the actual metric values.

    // line numbers in original file (= vector numbers) specifying which metric.
    int vector_num[NUM_WAY_MAX];
    vector_num[0] = atoi(argv_[1]);
    vector_num[1] = is_czek ? atoi(argv_[2]) : atoi(argv_[3]);
    vector_num[2] = num_way==2 ? 0 :
                    is_czek ? atoi(argv_[3]) : atoi(argv_[5]);
    assert(vector_num[0] >= 0 && vector_num[1] >= 0 && vector_num[2] >= 0);

    // table enrty index for this metric.
    int bit_num[NUM_WAY_MAX];
    bit_num[0] = is_czek ? 0 : atoi(argv_[2]);
    bit_num[1] = is_czek ? 0 : atoi(argv_[4]);
    bit_num[2] = num_way==2 ? 0 :
                 is_czek ? 0 : atoi(argv_[6]);

    int result = 0;
    if (is_czek && is_single)
      result = calculate_metric_elt<float, float>(metric_type_prec, num_way,
        vector_num, bit_num, snptxtfile, lifile);
    else if (is_czek) // && is_double
      result = calculate_metric_elt<double, double>(metric_type_prec, num_way,
        vector_num, bit_num, snptxtfile, lifile);
    else
      result = calculate_metric_elt<int, int>(metric_type_prec, num_way,
        vector_num, bit_num, snptxtfile, lifile);

    if (result != 0) {
      fprintf(stderr, "Error processing metric value.\n");
      return 1;
    }

    // Prepare for next line

    lineptr = 0;
    argc_ = 1;
  }

  // Finalize.

  fclose(snptxtfile);
  fclose(lifile);

  return 0;
}

//-----------------------------------------------------------------------------

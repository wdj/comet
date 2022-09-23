//-----------------------------------------------------------------------------
// "Manually" calculate metric values directly from original SNP (tped) file.
// This uses The CoMet output metrics file ONLY to get the required metrics
// coords, not the actual metrics values.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

enum {NUM_WAY_MAX = 3};

#define ASSERT(condition, command) { if (!(condition)) {command; return 1;} }

enum{MAX_LABEL_LEN = 255};

//-----------------------------------------------------------------------------

bool is_delim(int c) {
  return '\t' == c || ' ' == c;
}

//-----------------------------------------------------------------------------

template<class Field_t>
int calculate_metric_elt(char* metric_type, int num_way,
  int vector_num[NUM_WAY_MAX], int bit_num[NUM_WAY_MAX],
  FILE* snptxtfile, FILE* lifile) {

  // Initializations.

  const bool is_ccc = strcmp(metric_type, "ccc") == 0;
  const bool is_duo = strcmp(metric_type, "duo") == 0;
  const bool is_czek = strcmp(metric_type, "czekanowski") == 0;

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
  unsigned char* line = (unsigned char*)malloc(linesize);;

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

    // Loop to read the line from snptxtfile - loop up to newline.

    const size_t num_read_line = getline((char**)&line, &linesize, snptxtfile) - 1;
    ASSERT(line && num_read_line <= linesize,
      fprintf(stderr, "Error: memory %zu %zu\n", num_read_line, linesize));

    // Allocate elts[way_num] since now we have an upper bound on size.
    elts[way_num] = (Field_t*)malloc(num_read_line *
      num_allele_labels_per_field * sizeof(Field_t));

    for (size_t i=0, index=0, num_delim=0; i<num_read_line; ++i) {

      c = line[i];

      // Skip tab or space (these are separators).

      if (is_delim(c)) {
        num_delim++;
        continue;
      }

      // Get line label if at column 2.

      if (num_delim == 1) {
        // Append character to label.
        line_label[way_num][index++] = c; //check
        continue;
      }

      // If finished with frontmatter then reset char index.

      if (num_delim == num_frontmatter_fields - 1)
        index = 0;

      // Finished processing frontmatter.

      if (num_delim < num_frontmatter_fields)
        continue;

      // Store allele bit.

      elts[way_num][index++] = c;

      // If first cycle of num_way loop, then count num_field.

      if (0 == way_num && index % num_allele_labels_per_field == 0)
        num_field++;

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

    } // for i

  } // for way_num

  free(line);

  // Now we have all 2 (or 3) vectors and their associated data.

  // We will account for sparsity of the data







  // First get sum_i's and count_i's

  int count1[NUM_WAY_MAX];
  int sum1[NUM_WAY_MAX];

  for (int way_num=0; way_num<num_way; ++way_num) {
    count1[way_num] = 0;
    sum1[way_num] = 0;
    if (is_duo) {

      #pragma omp parallel for reduction(+:sum1[way_num]) reduction(+:count1[way_num])
      for (int f=0; f<num_field; ++f) {
        const int e0 = elts[way_num][f];

        if (e0 == '0')
          continue;

        count1[way_num] += 1;

        const int rho = (e0 == allele_labels[way_num][bit_num[way_num]]);
        sum1[way_num] += rho;
      } // for f

    } else { // if (!is_duo)

      for (int f=0; f<num_field; ++f) {
        const int e0 = elts[way_num][2*f];
        const int e1 = elts[way_num][2*f+1];

        if (e0 == '0')
          continue;

        count1[way_num] += 1;

        // Calculate row and add to sum - see paper

        const int rho = (e0 == allele_labels[way_num][bit_num[way_num]]) +
                        (e1 == allele_labels[way_num][bit_num[way_num]]);
        sum1[way_num] += rho;
      } // for f

    } // if (is_duo)

  } // for way_num

  // Now get sum_{ij}'s (or sum_{ijk}'s if 3-way)

  int countijk = 0;
  int sumijk = 0;

  if (is_duo) {

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

  } else { // if (!is_duo)

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

  } // if (is_duo)

  // substitute into formula

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

  const double value = 2 == num_way ?
    multiplier * fijk * (1 - param * f1[0]) *  (1 - param * f1[1]) :
    multiplier * fijk * (1 - param * f1[0]) *  (1 - param * f1[1])
                      * (1 - param * f1[2]);







  // Sort lines to output each result in a uniform order

  int perm[3];

  if (num_way == 2) {
    if (vector_num[0] > vector_num[1]) {
      perm[0] = 0;
      perm[1] = 1;
    } else {
      perm[0] = 1;
      perm[1] = 0;
    }
  } else {
    if (vector_num[0] > vector_num[1] && vector_num[1] > vector_num[2]) {
      perm[0] = 0;
      perm[1] = 1;
      perm[2] = 2;
    } else if (vector_num[0] > vector_num[2] && vector_num[2] > vector_num[1]) {
      perm[0] = 0;
      perm[1] = 2;
      perm[2] = 1;
    } else if (vector_num[1] > vector_num[0] && vector_num[0] > vector_num[2]) {
      perm[0] = 1;
      perm[1] = 0;
      perm[2] = 2;
    } else if (vector_num[1] > vector_num[2] && vector_num[2] > vector_num[0]) {
      perm[0] = 1;
      perm[1] = 2;
      perm[2] = 0;
    } else if (vector_num[2] > vector_num[0] && vector_num[0] > vector_num[1]) {
      perm[0] = 2;
      perm[1] = 0;
      perm[2] = 1;
    } else if (vector_num[2] > vector_num[1] && vector_num[1] > vector_num[0]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 0;
    }
  } // if num_way

  // Do output to stdout

  for (int way_num=0; way_num<num_way; ++way_num) {
    int iperm = perm[way_num];
    printf(0 != way_num ? " " : "");
    printf("%i %i", vector_num[iperm], bit_num[iperm]);
  } // for way_num

  for (int way_num=0; way_num<num_way; ++way_num) {
    int iperm = perm[way_num];
    printf(" %s_%c", line_label[iperm], allele_labels[iperm][bit_num[iperm]]);
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

// TODO metric_type_prec

  if (argc < 4) {
    printf("validate: create metrics data for validation of calculations\n");
    printf("Usage: validate "
           "<metric_type> <num_way> <snptxtfile> <line_index_file>\n");
    printf("Here stdin has the (ascii) metric entries, one per line.\n");
    return 0;
  }

  // Parse arguments.

  int argnum = 1;

  if (strcmp(argv[argnum], "ccc") != 0 &&
      strcmp(argv[argnum], "duo") != 0 &&
      strcmp(argv[argnum], "czekanowski") != 0) {
    fprintf(stderr, "Error: invalid metric_type. %s\n", argv[argnum]);
    return 1;
  }
  char* metric_type = argv[argnum++];
  const bool is_czek = strcmp(metric_type, "czekanowski") == 0;

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
  argnum++;

  FILE* lifile = fopen(argv[argnum], "rb");
  if (!lifile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

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





    int result = calculate_metric_elt<int>(metric_type, num_way, vector_num, bit_num,
                                           snptxtfile, lifile);

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

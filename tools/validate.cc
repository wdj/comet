//-----------------------------------------------------------------------------
// "Manually" calculate a CCC or DUO metric value directly from the original
// SNP file.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

enum{MAX_LABEL_LEN = 27};

//-----------------------------------------------------------------------------

int process_line(int argc, char** argv, const bool is_duo,
                 FILE* snptxtfile, FILE* lifile) {

  // Process args.

  if (argc != 1+4 && argc != 1+6) {
    // Input format expectd from each line of stdin:
    printf("Usage: "
      "validate_ccc_duo lineno0 bitno0 lineno1 bitno1 [lineno2 bitno2]\n");
    printf("Line and bit numbers are 0-based\n");
    return 1;
  }

  // Infer from number of command line args whether 2-way or 3-way.

  const int num_way = argc == 1+4 ? 2 : 3;
  enum {NUM_WAY_MAX = 3};

  // lineno values refer to the vectors being compared from the original file.
  int lineno[NUM_WAY_MAX];
  lineno[0] = atoi(argv[1]);
  lineno[1] = atoi(argv[3]);
  lineno[2] = num_way == 3 ? atoi(argv[5]) : 0;
  assert(lineno[0] >= 0);
  assert(lineno[1] >= 0);
  assert(lineno[2] >= 0);

  // bitno values refer to the table entry of the metric for those vectors.

  int bitno[NUM_WAY_MAX];
  bitno[0] = atoi(argv[2]);
  bitno[1] = atoi(argv[4]);
  bitno[2] = num_way == 3 ? atoi(argv[6]) : 0;
  assert(bitno[0] >= 0 && bitno[0] < 2);
  assert(bitno[1] >= 0 && bitno[1] < 2);
  assert(bitno[2] >= 0 && bitno[2] < 2);

  //CHECK
  // Allele labels would be "A", "C", "T", "G", and if none then use "X".
  // NOTE: we will use the term "allele label" to denote a single one of these letters.

  const int allele_label_null = 'X';

  // Initialize line labels to null.

  const int line_label_len = MAX_LABEL_LEN;
  unsigned char line_label[NUM_WAY_MAX][line_label_len+1];
  for (int i=0; i<num_way; ++i)
    for (int j=0; j<line_label_len+1; ++j)
      line_label[i][j] = 0;

  // Access the snp txt file.

  int fseek_success = fseek(snptxtfile, 0, SEEK_SET);
  if (0 != fseek_success) {
    fprintf(stderr, "Error: error reading SNP data file (0).\n");
    return 1;
  }

  // Allele label is represented by a character, stored as an int.

  typedef int AlleleLabel_t;

  AlleleLabel_t c = 0;

  // Read first line of snptxtfile to get upper bound on num fields per line.

  int num_read = 0;
  while ((c = fgetc(snptxtfile)) != '\n') {
    if (c == EOF) {
      fprintf(stderr, "Error snp file has no newlines.");
      return 1;
    }
    num_read++;
  }

  const int num_allele_labels_per_field = is_duo ? 1 : 2;

  // Round up. Note this is an overestimate because counting white space etc.
  const int num_field_max = (num_read + num_allele_labels_per_field - 1) / num_allele_labels_per_field;

  // Allocate storage for allele labels and allele bits for each needed SNP.

  enum {NUM_ALLELE_LABELS_PER_FIELD_MAX = 2};

  AlleleLabel_t allele_labels[NUM_WAY_MAX][NUM_ALLELE_LABELS_PER_FIELD_MAX];
  AlleleLabel_t* elts[NUM_WAY_MAX];
  for (int i=0; i<num_way; ++i) {
    allele_labels[i][0] = allele_label_null;
    allele_labels[i][1] = allele_label_null;
    elts[i] = (AlleleLabel_t*)malloc(num_field_max * NUM_ALLELE_LABELS_PER_FIELD_MAX * sizeof(AlleleLabel_t));
  }

  int num_field = 0;

  // Loop over num_way to input the required vectors.

  for (int way_num=0; way_num<num_way; ++way_num) {

    // Get location in the line index file needed to get the line index.

    fseek_success = fseek(lifile, lineno[way_num]*sizeof(size_t), SEEK_SET);
    if (0 != fseek_success) {
      fprintf(stderr, "Error: error reading line_indices file (1).\n");
      return 1;
    }

    // Get the index to the required line in the snptxtfile.

    size_t line_index = 0;
    num_read = fread(&line_index, sizeof(size_t), 1, lifile);
    if (1 != num_read) {
      fprintf(stderr, "Error: error reading line_indices file (2).\n");
      return 1;
    }

    // Set next-read pointer to start of needed line in snptxtfile.

    fseek_success = fseek(snptxtfile, line_index, SEEK_SET);
    if (0 != fseek_success) {
      fprintf(stderr, "Error: error reading SNP data file (1).\n");
      return 1;
    }

    const int num_frontmatter_fields = 4;
    int num_delim = 0;
    int index = 0;

    // Loop to read the line from snptxtfile - loop up to newline.

    while ((c = fgetc(snptxtfile)) != '\n') {

      // Skip tab or space (these are separators).

      if (c == '\t' || c == ' ') {
        num_delim++;
        continue;
      }

      // Get line label

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

      if (way_num == 0)
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

    } // while c

  } // for way_num

  num_field /= num_allele_labels_per_field;

  // Now we have all 2 (or 3) vectors and their associated data.

  // We will account for sparsity of the data

  // First get sum_i's and count_i's

  int count1[NUM_WAY_MAX];
  int sum1[NUM_WAY_MAX];

  for (int way_num=0; way_num<num_way; ++way_num) {
    count1[way_num] = 0;
    sum1[way_num] = 0;
    if (is_duo) {

      for (int f=0; f<num_field; ++f) {
        const int e0 = elts[way_num][f];

        if (e0 == '0')
          continue;

        count1[way_num] += 1;

        const int rho = (e0 == allele_labels[way_num][bitno[way_num]]);
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

        const int rho = (e0 == allele_labels[way_num][bitno[way_num]]) +
                        (e1 == allele_labels[way_num][bitno[way_num]]);
        sum1[way_num] += rho;
      } // for f

    } // if (is_duo)

  } // for way_num

  // Now get sum_{ij}'s (or sum_{ijk}'s if 3-way)

  int countijk = 0;
  int sumijk = 0;

  if (is_duo) {

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

      const int rho0 = (e00 == allele_labels[0][bitno[0]]);
      const int rho1 = (e10 == allele_labels[1][bitno[1]]);
      const int rho2 = num_way == 2 ? 1 :
                       (e20 == allele_labels[2][bitno[2]]);

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

      const int rho0 = (e00 == allele_labels[0][bitno[0]]) +
                       (e01 == allele_labels[0][bitno[0]]);
      const int rho1 = (e10 == allele_labels[1][bitno[1]]) +
                       (e11 == allele_labels[1][bitno[1]]);
      const int rho2 = num_way == 2 ? 1 :
                       (e20 == allele_labels[2][bitno[2]]) +
                       (e21 == allele_labels[2][bitno[2]]);

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
    if (lineno[0] > lineno[1]) {
      perm[0] = 0;
      perm[1] = 1;
    } else {
      perm[0] = 1;
      perm[1] = 0;
    }
  } else {
    if (lineno[0] > lineno[1] && lineno[1] > lineno[2]) {
      perm[0] = 0;
      perm[1] = 1;
      perm[2] = 2;
    } else if (lineno[0] > lineno[2] && lineno[2] > lineno[1]) {
      perm[0] = 0;
      perm[1] = 2;
      perm[2] = 1;
    } else if (lineno[1] > lineno[0] && lineno[0] > lineno[2]) {
      perm[0] = 1;
      perm[1] = 0;
      perm[2] = 2;
    } else if (lineno[1] > lineno[2] && lineno[2] > lineno[0]) {
      perm[0] = 1;
      perm[1] = 2;
      perm[2] = 0;
    } else if (lineno[2] > lineno[0] && lineno[0] > lineno[1]) {
      perm[0] = 2;
      perm[1] = 0;
      perm[2] = 1;
    } else if (lineno[2] > lineno[1] && lineno[1] > lineno[0]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 0;
    }
  } // if num_way

  // Do output to stdout

  for (int i=0; i<num_way; ++i) {
    int iperm = perm[i];
    printf(0 != i ? " " : "");
    printf("%i %i", lineno[iperm], bitno[iperm]);
  } // for i

  for (int i=0; i<num_way; ++i) {
    int iperm = perm[i];
    printf(" %s_%c", line_label[iperm], allele_labels[iperm][bitno[iperm]]);
  } // for i

  printf(" %f\n", value);

  for (int i=0; i<num_way; ++i) {
    free(elts[i]);
  }
  return 0;
}

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  // Help message.

  if (argc < 4) {
    printf("validate_ccc_duo: create validation data for calculations\n");
    printf("Usage: validate_ccc_duo <metric_type> <num_way> <snptxtfile> <line_index_file>\n");
    printf("Here stdin is composed of metric entries, one per line.\n");
    return 0;
  }

  // Parse arguments.

  if (strcmp(argv[1], "ccc") != 0 && strcmp(argv[1], "duo") != 0) {
    fprintf(stderr, "Error: invalid metric_type\n");
    return 1;
  }

  if (strcmp(argv[2], "2") != 0 && strcmp(argv[2], "3") != 0) {
    fprintf(stderr, "Error: invalid num_way\n");
    return 1;
  }

  const bool is_duo = strcmp(argv[1], "duo") == 0;
  const int num_way = atoi(argv[2]);
  char* snptxtfilepath = argv[3];
  char* lifilepath = argv[4];

  // Access input files.

  FILE* snptxtfile = fopen(snptxtfilepath, "r");
  if (!snptxtfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", snptxtfilepath);
    return 1;
  }
//  int fseek_success = fseek(snptxtfile, 0, SEEK_SET);
//  if (0 != fseek_success) {
//    fprintf(stderr, "xxxError: error reading SNP data file (0).\n");
//    return 1;
//  }

  FILE* lifile = fopen(lifilepath, "rb");
  if (!lifile) {
    fprintf(stderr, "Error: unable to open file. %s\n", lifilepath);
    return 1;
  }

#if 0
for (size_t line_num=0; line_num<10 ; ++line_num) {

    int fseek_success = fseek(lifile, line_num*sizeof(size_t), SEEK_SET);
    if (0 != fseek_success) {
      fprintf(stderr, "Error: error reading line_indices file (1).\n");
      return 1;
    }

    // Get the index to the required line in the snptxtfile.

    size_t line_index = 0;
    int num_read = fread(&line_index, sizeof(size_t), 1, lifile);
    if (1 != num_read) {
      fprintf(stderr, "Error: error reading line_indices file (xx).\n");
      return 1;
    }

printf("%zu %zu\n", line_num, line_index);
}
#endif

  // Prepare to read from stdin.

  enum{LINE_LEN_MAX = 4096};
  unsigned char line[LINE_LEN_MAX];

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

    // Create an (argc, argv) arg list consisting of the tokens.

    line[lineptr] = 0;
    argv_[0] = (char*)&line[lineptr];
    for (int i=0; i<lineptr; ++i) {
      if (line[i] == ' ') {
        line[i] = 0;
      }
    }
    for (int i=0; i<lineptr; ++i) {
      if (line[i] != 0 && (i == 0 || line[i-1] == 0)) {
          argv_[argc_++] = (char*)&line[i];
      }
    }
    argc_ = argc_ < 2*num_way+1 ? argc_ : 2*num_way+1;

    // Process this line.

    // Only use the first several tokens, specifying line and bit numbers;
    // discard the actual metric values.

    int result = process_line(argc_, argv_, is_duo, snptxtfile, lifile);
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

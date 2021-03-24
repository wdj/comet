//-----------------------------------------------------------------------------
// Postprocess a single output file from CoMet - convert binary to text.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//enum {MAX_VECTOR_LABEL_LEN = 27};
enum {MAX_VECTOR_LABEL_LEN = 255};

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  if (argc < 6) {
    printf("postprocess_ccc_duo: convert metrics file from binary to text\n");
    printf("Usage: postprocess <num_way> <allele_label_file>"
           " <vector_label_file> <metrics_bin_file> <metrics_txt_file>\n");
    return 0;
  }

  int argnum = 1;

  //if (strcmp(argv[argnum], "ccc") != 0 && strcmp(argv[argnum], "duo") != 0) {
  //  fprintf(stderr, "Error: invalid metric_type\n");
  //  return 1;
  //}
  //const bool is_duo = strcmp(argv[argnum++], "duo") == 0;

  if (!( ('2' == argv[argnum][0] || '3' == argv[argnum][0]) &&
         0 == argv[argnum][1] )) {
    fprintf(stderr, "Error: invalid num_way\n");
    return 1;
  }
  const int num_way = atoi(argv[argnum++]);

  FILE* allele_label_file = fopen(argv[argnum], "r");
  if (!allele_label_file) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  FILE* vector_label_file = fopen(argv[argnum], "r");
  if (!vector_label_file) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  char* metricsbinfilename = argv[argnum];
  FILE* metricsbinfile = fopen(argv[argnum], "rb");
  if (!metricsbinfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  FILE* metricstxtfile = fopen(argv[argnum], "w");
  if (!metricstxtfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  // Initializations

  enum {NUM_WAY_MAX = 3};

  int num_read = 0;
  int fseek_success = 0;
  float value = 0;
  int coord[NUM_WAY_MAX];
  const int no_allele = 'X'; // Marker for case only one allele label in a line

  // Determine vector label length by examining first vector label.

  int i = 0;
  for (i=0; i<MAX_VECTOR_LABEL_LEN; ++i) {

    fseek_success = fseek(vector_label_file, i, SEEK_SET);
    if (0 != fseek_success) {
      fprintf(stderr, "Error: error reading file. 0.1\n");
      return 1;
    }

    unsigned char c = 0;
    const int num_to_read = 1;

    num_read = fread(&c, sizeof(unsigned char), num_to_read, vector_label_file);
    if (num_to_read != num_read) {
      fprintf(stderr, "Error: error reading file. 0.2\n");
      return 1;
    }

    if (c == '\n')
      break;

  } // for i

  if (MAX_VECTOR_LABEL_LEN == i) {
    fprintf(stderr, "Error in vector labels file.\n");
    return 1;
  }

  const int vector_label_len = i;
  const int allele_label_len = 2;

  // Loop over result values in specified file
  // while statement reads first coord.

  while ((num_read = fread(&coord[0], sizeof(int), 1, metricsbinfile)) == 1) {

    // Get (rest of) coordinates and metric value

    num_read = fread(&coord[1], sizeof(int), 1, metricsbinfile);
    if (1 != num_read) {
      fprintf(stderr, "Error: error reading file. 1\n");
      return 1;
    }

    if (3 == num_way) {
      num_read = fread(&coord[2], sizeof(int), 1, metricsbinfile);
      if (num_read != 1) {
        fprintf(stderr, "Error: error reading file. 2\n");
        return 1;
      }
    }

    num_read = fread(&value, sizeof(float), 1, metricsbinfile);
    if (1 != num_read) {
      fprintf(stderr, "Error: error reading file. 3\n");
      return 1;
    }

    int lineno[NUM_WAY_MAX];
    int bitno[NUM_WAY_MAX];
    unsigned char allele_label[NUM_WAY_MAX][allele_label_len];
    unsigned char vector_label[NUM_WAY_MAX][vector_label_len+1];

    // Loop over coordinates

    for (int way_num=0; way_num<num_way; ++way_num) {

      // Decode vector number (line number in file) and whether lo or hi bit

      lineno[way_num] = coord[way_num] / 2;
      bitno[way_num] = coord[way_num] % 2;

      // Get allele labels

      fseek_success = fseek(allele_label_file,
        (allele_label_len+1)*lineno[way_num], SEEK_SET);
      if (0 != fseek_success) {
        fprintf(stderr, "Error: error reading file. 4\n");
        return 1;
      }

      num_read = fread(&allele_label[way_num], sizeof(unsigned char),
        allele_label_len, allele_label_file);
      if (allele_label_len != num_read) {
        fprintf(stderr, "Error: error reading file. 5 %s\n",
                metricsbinfilename);
        return 1;
      }

      // If repeated then disregard upper value, replace with "X"
      if (allele_label[way_num][0] == allele_label[way_num][1]) {
        allele_label[way_num][1] = no_allele;
      }

      // Get vector label

      fseek_success = fseek(vector_label_file,
        (vector_label_len+1)*lineno[way_num], SEEK_SET);
      if (0 != fseek_success) {
        fprintf(stderr, "Error: error reading file. 6\n");
        return 1;
      }

      num_read = fread(&vector_label[way_num], sizeof(unsigned char),
                       vector_label_len, vector_label_file);
      if (vector_label_len != num_read) {
        fprintf(stderr, "Error: error reading file. 7\n");
        return 1;
      }

      // Remove end padding
      for (int j=0; j<vector_label_len; ++j) {
        if (' ' == (int)vector_label[way_num][j]) {
          vector_label[way_num][j] = 0;
        }
      }
      vector_label[way_num][vector_label_len] = 0;

    } // for way_num

    // Permute labels to output each result with a uniform order of labels
    // By convention let line numbers bein increasing, e.g. "0 1" not "1 0"

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

    // Output line and bit numbers

    for (int i=0; i<num_way; ++i) {
      int iperm = perm[i];
      fprintf(metricstxtfile, 0 != i ? " " : "");
      fprintf(metricstxtfile, "%i %i", lineno[iperm], bitno[iperm]);
    } // for i

    // Output vector label

    for (int way_num=0; way_num<num_way; ++way_num) {
      int iperm = perm[way_num];
      fprintf(metricstxtfile, " %s_%c", vector_label[iperm], allele_label[iperm][bitno[iperm]]);
    } // for i

    // Output value

    fprintf(metricstxtfile, " %f\n", value);

  } // while

  fclose(metricsbinfile);
  fclose(allele_label_file);
  fclose(vector_label_file);
  fclose(metricstxtfile);

  return 0;
}

//-----------------------------------------------------------------------------

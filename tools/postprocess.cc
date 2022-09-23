//-----------------------------------------------------------------------------
// Postprocess a single output file from CoMet - convert binary to text.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  // Parse arguments.

  if (argc < 6) {
    printf("postprocess: convert metrics file from binary to text\n");
    printf("Usage: postprocess [<metric_type>] <num_way> <allele_label_file>"
           " <line_label_file> <metrics_bin_file> <metrics_txt_file>\n");
    return 0;
  }

  int argnum = 1;

  bool is_czek;
  if (7 == argc) {
    if (strcmp(argv[argnum], "ccc") != 0 &&
        strcmp(argv[argnum], "duo") != 0 &&
        strcmp(argv[argnum], "czekanowski") != 0) {
      fprintf(stderr, "Error: invalid metric_type. %s\n", argv[argnum]);
      return 1;
    }
    is_czek = strcmp(argv[argnum++], "czekanowski") == 0;
  }

  if (!( ('2' == argv[argnum][0] || '3' == argv[argnum][0]) &&
         0 == argv[argnum][1] )) {
    fprintf(stderr, "Error: invalid num_way\n");
    return 1;
  }
  const int num_way = atoi(argv[argnum++]);

  FILE* allele_label_file = is_czek ? NULL : fopen(argv[argnum], "r");
  if (!is_czek && !allele_label_file) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  FILE* line_label_file = fopen(argv[argnum], "r");
  if (!line_label_file) {
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
  int coord[NUM_WAY_MAX];
  const int no_allele = 'X'; // Marker for case only one allele label in a line

  // Determine line label length by examining first line label.
  // NOTE: this is coordinated with line_labels.sh - uniform label lengths.

  int i = 0;
  for (i=0 ;; ++i) {

    const int fseek_success = fseek(line_label_file, i, SEEK_SET);
    if (0 != fseek_success) {
      fprintf(stderr, "Error: error reading file. 0.1\n");
      return 1;
    }

    unsigned char c = 0;
    const int num_to_read = 1;

    int num_read = fread(&c, sizeof(unsigned char), num_to_read, line_label_file);
    if (num_to_read != num_read) {
      fprintf(stderr, "Error: error reading file. 0.2\n");
      return 1;
    }

    if (c == '\n')
      break;

  } // for i

  const int line_label_len = i;
  const int allele_label_len = 2;

  // LOOP over result values in specified file.

  // while statement reads first coord of metric.

  while ((num_read = fread(&coord[0], sizeof(int), 1, metricsbinfile)) == 1) {

    // Get (rest of the) coordinates and metric value

    for (int i=1; i<num_way; ++i) {
      const int num_read = fread(&coord[i], sizeof(int), 1, metricsbinfile);
      if (1 != num_read) {
        fprintf(stderr, "Error: error reading metric coord\n");
        return 1;
      }
    }

//    num_read = fread(&coord[1], sizeof(int), 1, metricsbinfile);
//    if (1 != num_read) {
//      fprintf(stderr, "Error: error reading file. 1\n");
//      return 1;
//    }
//
//    if (3 == num_way) {
//      num_read = fread(&coord[2], sizeof(int), 1, metricsbinfile);
//      if (num_read != 1) {
//        fprintf(stderr, "Error: error reading file. 2\n");
//        return 1;
//      }
//    }

    float value = 0;
    int num_read = fread(&value, sizeof(float), 1, metricsbinfile);
    if (1 != num_read) {
      fprintf(stderr, "Error: error reading file. 3\n");
      return 1;
    }

    int lineno[NUM_WAY_MAX];
    int bitno[NUM_WAY_MAX];
    unsigned char allele_label[NUM_WAY_MAX][allele_label_len];
    unsigned char line_label[NUM_WAY_MAX][line_label_len+1];

    // LOOP over coordinates.

    for (int way_num=0; way_num<num_way; ++way_num) {

      // Decode vector number (line number in file) and whether lo or hi bit.

      const int num_bitno = is_czek ? 1 : 2;

      lineno[way_num] = coord[way_num] / num_bitno;
      bitno[way_num] = coord[way_num] % num_bitno;

      // Get allele labels if needed.

      allele_label[way_num][0] == no_allele;
      allele_label[way_num][1] == no_allele;
      if (allele_label_file) {
        const int fseek_success = fseek(allele_label_file,
          (allele_label_len+1)*lineno[way_num], SEEK_SET);
        if (0 != fseek_success) {
          fprintf(stderr, "Error: error reading file. 4\n");
          return 1;
        }

        const int num_read = fread(&allele_label[way_num],
          sizeof(unsigned char), allele_label_len, allele_label_file);
        if (allele_label_len != num_read) {
          fprintf(stderr, "Error: error reading file. 5 %s  %i\n",
                  metricsbinfilename, lineno[way_num]);
          return 1;
        }
      }

      // If repeated then disregard upper value, replace with "X"
      if (allele_label[way_num][0] == allele_label[way_num][1])
        allele_label[way_num][1] = no_allele;

      // Get line label.

      fseek_success = fseek(line_label_file,
        (line_label_len+1)*lineno[way_num], SEEK_SET);
      if (0 != fseek_success) {
        fprintf(stderr, "Error: error reading file. 6\n");
        return 1;
      }

      num_read = fread(&line_label[way_num], sizeof(unsigned char),
                       line_label_len, line_label_file);
      if (line_label_len != num_read) {
        fprintf(stderr, "Error: error reading file. 7\n");
        return 1;
      }

      // Remove line label end padding that was added to help with indexing.
      for (int j=0; j<line_label_len; ++j) {
        if (' ' == (int)line_label[way_num][j]) {
          line_label[way_num][j] = 0;
        }
      }
      line_label[way_num][line_label_len] = 0;

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

    for (int way_num=0; way_num<num_way; ++way_num) {
      int iperm = perm[way_num];
      fprintf(metricstxtfile, 0 != way_num ? " " : "");
      fprintf(metricstxtfile, "%i", lineno[iperm]);
      if (!is_czek)
        fprintf(metricstxtfile, " %i", bitno[iperm]);
    } // for way_num

    // Output line label

    for (int way_num=0; way_num<num_way; ++way_num) {
      int iperm = perm[way_num];
      fprintf(metricstxtfile, " %s", line_label[iperm]);
      if (!is_czek)
        fprintf(metricstxtfile, "_%c", allele_label[iperm][bitno[iperm]]);
    } // for i

    // Output value

    fprintf(metricstxtfile, " %f\n", value);

  } // while

  fclose(metricsbinfile);
  if (!is_czek)
    fclose(allele_label_file);
  fclose(line_label_file);
  fclose(metricstxtfile);

  return 0;
}

//-----------------------------------------------------------------------------

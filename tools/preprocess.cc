//-----------------------------------------------------------------------------
// Convert a SNP file in tped format to a binary file.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  // Parse arguments.

  if (argc != 4) {
    printf("preprocess: convert text SNP file to a packed binary SNP file.\n");
    printf("Usage: preprocess <metric_type> <snp_text_file> <snp_bin_file>\n");
    return 0;
  }

  int argnum = 1;

  if (strcmp(argv[argnum], "ccc") != 0 && strcmp(argv[argnum], "duo") != 0) {
    fprintf(stderr, "Error: invalid metric_type\n");
    return 1;
  }
  const bool is_duo = strcmp(argv[argnum], "duo") == 0;
  argnum++;

  FILE* snptxtfile = fopen(argv[argnum], "r");
  if (!snptxtfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  FILE* snpbinfile = fopen(argv[argnum], "wb");
  if (!snpbinfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  // Initializations.

  const int max_line_len = 1000000;
  int* line = (int*)malloc(max_line_len * sizeof(*line));

  const int num_frontmatter_fields = 4;

  int c = 0;
  int cprev = 0;
  int line_len = 0;
  int clo = 0;
  int chi = 0;
  int col = 0;

  // Loop to input chars from input file

  while ((c = fgetc(snptxtfile)) != EOF) {

    line[line_len] = c;
    line_len++;
    if (line_len > max_line_len) {
      fprintf(stderr, "Error: input line too long\n");
      return 1;
    }

    if (c != '\n') {
      continue;
    }

    // Reached end of a line, so process it

    line_len--;

    //----------
    // PASS 1: loop over tokens to get the allele labels
    //----------

    int i, num_delims;
    for (i=0, num_delims=0; i<line_len; ++i) {

      c = line[i];
      if (c == '\t' || c == ' ') {
        // Found delimiter, go to next token
        num_delims++;
        continue;
      }
      if (num_delims < num_frontmatter_fields) {
        // Skip first four tokens - pick these up with another command
        continue;
      }

      // Get token number
      col = num_delims - num_frontmatter_fields;

      if (col % 2 == 1 && ! is_duo) {
        if ((c=='0' && cprev!='0') ||
            (c!='0' && cprev=='0')) {
          fprintf(stderr, "Error: token pair must be both zero or both nonzero\n");
          return 1;
        }
      }

      // Record values of the tokens encountered

      if (c == '0') {
      } else if (clo == 0) {
        clo = c;
      } else if (clo < c && chi == 0) {
        chi = c;
      } else if (clo > c && chi == 0) {
        chi = clo;
        clo = c;
      } else if (c != clo && c != chi) {
        fprintf(stderr, "Error: malformed input snp file.\n");
        return 1;
      }

      cprev = c;
    } // for i - first pass loop

    //----------

    if (0 == chi)
      chi = clo;

    if (col % 2 == 0 && ! is_duo) {
      fprintf(stderr, "Error: line has invalid number of tokens\n");
      return 1;
    }

    int out_buf = 0;
    int num_buf = 0;

    //----------
    // PASS 2: loop to output results
    //----------

    for (i=0, num_delims=0; i<line_len; ++i) {

      c = line[i];
      if (c == '\t' || c == ' ') {
        // Found delimiter, go to next token
        num_delims++;
        continue;
      }
      if (num_delims < num_frontmatter_fields) {
        // Skip first four tokens - pick these up with another command
        continue;
      }

      // Get token number
      col = num_delims - num_frontmatter_fields;

      if (col % 2 == 1 && ! is_duo) {
        cprev = c;
        continue;
      }

      const int c0 = cprev;
      const int c1 = c;

      // Next map the token pair to a pair of bits needed by CoMet.

      // Note the alphabetically first allele is mapped to the zero bit

      // NOTE: order of lines matters below.
      int sn = is_duo ? (
            c1 == '0' ? 2 * (1) + 1 * (0)  //  0  ->  1 0
          : c1 == clo ? 2 * (0) + 1 * (0)  //  A  ->  0 0
          : c1 == chi ? 2 * (0) + 1 * (1)  //  B  ->  0 1
          : -1
        ) : (
            c0 == '0' && c1 == '0' ? 2 * (1) + 1 * (0)  //  00  ->  1 0
          : c0 == clo && c1 == clo ? 2 * (0) + 1 * (0)  //  AA  ->  0 0
          : c0 == clo && c1 == chi ? 2 * (0) + 1 * (1)  //  AB  ->  0 1
          : c0 == chi && c1 == clo ? 2 * (0) + 1 * (1)  //  BA  ->  0 1
          : c0 == chi && c1 == chi ? 2 * (1) + 1 * (1)  //  BB  ->  1 1
          : -1
        );

      if (sn == -1) {
        fprintf(stderr, "Error: unknown error encountered mapping tokens\n");
        return 1;
      }

      // Insert the 2 bits into the output buffer

      out_buf = (out_buf << 2) | sn;
      num_buf++;

      // Output buffer has 8 bits, so flush

      if (num_buf == 4) {
        fputc(out_buf, snpbinfile);
        out_buf = 0;
        num_buf = 0;
      }

    } // for i - second pass loop

    //----------

    // Final flush of buffer

    if (num_buf != 0) {
      out_buf = out_buf << (2 * (4 - num_buf));
      fputc(out_buf, snpbinfile);
      out_buf = 0;
      num_buf = 0;
    }

    // Process next line

    line_len = 0;
    clo = 0;
    chi = 0;

  } // while c

  // Finish

  free(line);

  fclose(snptxtfile);
  fclose(snpbinfile);

  return 0;
}

//-----------------------------------------------------------------------------

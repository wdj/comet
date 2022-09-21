//-----------------------------------------------------------------------------
// Validate the result of running the preprocess command.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {

  if (argc != 5) {
    printf("preprocess_validate: check binary SNP file against text SNP file.\n");
    printf("Usage: preprocess_validate <metric_type> <snp_bin_file> <snp_text_file> <allele_label_file>\n");
    return 0;
  }

  int argnum = 1;

  char* metric_type = argv[argnum];
  
  if (strcmp(metric_type, "ccc") != 0 && strcmp(metric_type, "duo") != 0) {
    fprintf(stderr, "Error: invalid metric_type\n");
    return 1; 
  }      
  const bool is_duo = strcmp(argv[argnum], "duo") == 0;
  argnum++;

  FILE* snpbinfile = fopen(argv[argnum], "rb");
  if (!snpbinfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  } 
  argnum++;

  FILE* snptxtfile = fopen(argv[argnum], "r");
  if (!snptxtfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  } 
  argnum++;
  
  FILE* alfile = fopen(argv[argnum], "r");
  if (!alfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  } 
  argnum++;

  // Initializations.

  const int num_frontmatter_fields = 4;

  size_t num_checked = 0;
  size_t num_validated = 0;

  // Loop over lines

  for (int line_num=0; ; ++line_num) {

    // Read from allele label file to get info for this line

    int c = 0;
    if ((c = fgetc(alfile)) == EOF) {
      break; // Check complete
    }

    const int clo = c;
    const int chi = fgetc(alfile);
    int discard = fgetc(alfile); // skip newline

    // Skip first tokens in line of snptxtfile

    int num_delims = 0;
    while ((c = fgetc(snptxtfile)) != EOF) {
      if (c == '\t' || c == ' ') {
        num_delims++;
      }
      if (num_delims == num_frontmatter_fields) {
        break;
      }
    } // while

    // Loop over elements of this line

    int inbyte = 0;
    for (int elt_num=0; ; ++elt_num) {

      // Pick up true allele information from snp text file.

      const int c0true = fgetc(snptxtfile);
      discard = fgetc(snptxtfile); // skip delim
      const int c1true = is_duo ? '0' : fgetc(snptxtfile);
      if (! is_duo)
        discard = fgetc(snptxtfile); // skip delim
      const int is_end_of_line = discard == 10;

      // Extract 2-bit seminibble from snp bin file

      const int sn_num = elt_num % 4;

      if (sn_num == 0)
        inbyte = fgetc(snpbinfile); // Get new byte

      const int sn = (inbyte >> (6 - 2*sn_num) ) & 3;
//printf("%i %i\n", elt_num, sn);

      // Map

      int c0 = 0;
      int c1 = 0;

      if (is_duo) {
        c1 = '0';
        if (sn == 0) {
          c0 = clo;
        } else if (sn == 1) {
          c0 = chi;
        } else if (sn == 2) {
          c0 = '0';
        } else if (sn == 3) {
          c0 = 'X'; // should never happen.
        }
      } else { // ! is_duo
        if (sn == 0) {
          c0 = clo;
          c1 = clo;
        } else if (sn == 1) {
          c0 = clo;
          c1 = chi;
        } else if (sn == 2) {
          c0 = '0';
          c1 = '0';
        } else if (sn == 3) {
          c0 = chi;
          c1 = chi;
        }
      } // if (is_duo)

      // Check

      if ((c0==c0true && c1==c1true) || (c0==c1true && c1==c0true)) {
        num_validated++;
      } else {
        printf("Error: invalid value detected, line %i elt %i\n", line_num+1, elt_num+1);
        if (is_duo)
          printf("actual %c  expected %c\n", c0, c0true);
        else
          printf("actual %c %c  expected %c %c\n", c0, c1, c0true, c1true);
        return 1;
      }
      num_checked++;

      if (is_end_of_line) {
        break;
      }

    } // elt_num

  } // line_num

  printf("%s metric number of elements validated: %ul of %ul\n",
         metric_type, num_validated, num_checked);

  // Finish

  fclose(alfile);
  fclose(snpbinfile);
  fclose(snptxtfile);

  return 0;
}

//-----------------------------------------------------------------------------

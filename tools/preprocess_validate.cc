//-----------------------------------------------------------------------------
// Validate the result of running the preprocess command.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>

//-----------------------------------------------------------------------------

bool is_delim(int c) {
  return '\t' == c || ' ' == c;
}

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  // Parse arguments.

  if (!(argc == 5 || (argc == 4 && (
          strcmp(argv[1], "czekanowski_single") == 0 ||
          strcmp(argv[1], "czekanowski_double") == 0)))) {
    printf("preprocess_validate: check binary SNP (CoMet input) file "
           "against text SNP (tped) file.\n");
    printf("Usage: preprocess_validate <metric_type_prec> <snp_bin_file> "
           "<snp_text_file> [<allele_label_file>]\n");
    return 0;
  }

  int argnum = 1;
  char* metric_type_prec = argv[argnum];

  if (strcmp(metric_type_prec, "ccc") != 0 &&
      strcmp(metric_type_prec, "duo") != 0 &&
      strcmp(metric_type_prec, "czekanowski_single") != 0 &&
      strcmp(metric_type_prec, "czekanowski_double") != 0) {
    fprintf(stderr, "Error: invalid metric_type_prec\n");
    return 1;
  }
  const bool is_ccc = strcmp(metric_type_prec, "ccc") == 0;
  const bool is_duo = strcmp(metric_type_prec, "duo") == 0;
  const bool is_czek = strcmp(metric_type_prec, "czekanowski_single") == 0 ||
                       strcmp(metric_type_prec, "czekanowski_double") == 0;
  const bool is_single = strcmp(metric_type_prec, "czekanowski_single") == 0;

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
  FILE* alfile = !is_czek ? fopen(argv[argnum], "r") : NULL;
  if (!is_czek) {
    if (!alfile) {
      fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
      return 1;
    } 
  } 
  argnum++;

  // Initializations.

  const int num_frontmatter_fields = 4;

  size_t num_checked = 0;
  size_t num_validated = 0;

  enum {NULL_AL = -1};

  //----------
  // LOOP over lines.
  //----------

  for (int line_num=0; ; ++line_num) {

    int c = 0;
    int clo = NULL_AL;
    int chi = NULL_AL;

    // If appropriate, read from allele label file to get info for this line.

    if (!is_czek) {
      if ((c = fgetc(alfile)) != EOF) {
        clo = c;
        if ((c = fgetc(alfile)) != EOF) {
          chi = c;
          if ((c = fgetc(alfile)) != EOF) {
            // Skip newline.
          }
        }
      }
    }

    // Skip first tokens in line of snptxtfile

    int num_delims = 0;
    while ((c = fgetc(snptxtfile)) != EOF) {
      if (is_delim(c))
        ++num_delims;
      if (num_delims == num_frontmatter_fields)
        break;
    } // while
    if (EOF == c)
      break;

    enum {NEWLINE = 10};

    //----------
    // LOOP over elements of this line.
    //----------

    int input_buf = 0;
    for (int elt_num=0; ; ++elt_num) {

      //----------
      // 1. Handle ccc, duo case.
      //----------

      if (!is_czek) {

        // Pick up true element information from text file.
        // NOTE: assuming here that the tped file is syntactically correct.

        const int c0true = fgetc(snptxtfile);
        int discard = fgetc(snptxtfile); // skip delim
        const int c1true = is_duo ? '0' : fgetc(snptxtfile);
        if (!is_duo)
          discard = fgetc(snptxtfile); // skip delim
        const int is_end_of_line = discard == NEWLINE;

        // Extract element from bin file (2 bits).

        const int sn_num = elt_num % 4;

        if (sn_num == 0) {
          // Get new byte.
          if((input_buf = fgetc(snpbinfile)) == EOF) {
            fprintf(stderr, "Error: premature termination of bin file, "
                            "line %i elt %i\n", line_num+1, elt_num+1);
            return 1;
          }
        }

        //const int sn = (input_buf >> (6 - 2*sn_num) ) & 3;
        const int sn = (input_buf >> (2*sn_num) ) & 3;

//printf("%i %i %i %i %i\n", sn_num, sn, input_buf, elt_num, is_end_of_line);

        // Map snp value to c0/c1.

        int c0 = 0;
        int c1 = 0;

        if (is_duo) {
          c1 = '0';
          if (sn == 0)
            c0 = clo;
          else if (sn == 1)
            c0 = chi;
          else if (sn == 2)
            c0 = '0';
          else if (sn == 3)
            c0 = 'X'; // should never happen.
        } else { // is_ccc
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

        // Check.

        if ((c0==c0true && c1==c1true) || (c0==c1true && c1==c0true))
          num_validated++;
        else {
          fprintf(stderr, "Error: invalid value detected, line %i elt %i\n",
                  line_num+1, elt_num+1);
          if (is_duo)
            fprintf(stderr, "actual %c  expected %c\n", c0, c0true);
          else
            fprintf(stderr, "actual %c %c  expected %c %c\n", c0, c1,
                    c0true, c1true);
          return 1;
        }

        num_checked++;

        if (is_end_of_line)
          break;

      } else { // is_czek

      //----------
      // 2. Handle Czekanowski case.
      //----------

        // Pick up true element information from text file.
        // NOTE: assuming here that the tped file is syntactically correct.

        std::string token;

        while ((c = fgetc(snptxtfile)) != EOF && !is_delim(c) && c != NEWLINE)
          token += c;
        const int is_end_of_line = c == NEWLINE;

        // Extract element from bin file (4 or 8 bytes).

        void* value = 0;
        //assert(sizeof(void*)>=sizeof(float) && sizeof(void*)>=sizeof(double));

        const size_t num_read = fread(&value,
          is_single ? sizeof(float) : sizeof(double), 1, snpbinfile);

        if (1 != num_read) {
          fprintf(stderr, "Error: premature termination of bin file.\n");
          return 1;
        }

        // Check.

        if ( (is_single && *(float*)&value != std::stof(token)) || (
             !is_single && *(double*)&value != std::stod(token)) ) {
          fprintf(stderr, "Error: invalid value detected, line %i elt %i\n",
                  line_num+1, elt_num+1);
          return 1;
        } else
          num_validated++;

        num_checked++;

        if (is_end_of_line)
          break;

      } // !is_czek

      //----------
      // End handling of cases.
      //----------

    } // elt_num

  } // line_num

  printf("%s metric number of elements validated: %ul of %ul: %s\n",
         metric_type_prec, num_validated, num_checked,
         num_validated == num_checked ? "=== PASSED ===" : "");

  // Finish

  if (alfile)
    fclose(alfile);
  fclose(snpbinfile);
  fclose(snptxtfile);

  return 0;
}

//-----------------------------------------------------------------------------

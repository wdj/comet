//-----------------------------------------------------------------------------
// Convert a SNP file in tped format to a binary file.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <string>
#include <sstream>
#include <iostream>

//-----------------------------------------------------------------------------

class ElementWriter {

  // 4 seminibbles per byte.
  enum {BUF_SIZE = 4};

public:

  ElementWriter(char* const output_file_name)
    : output_file_(fopen(output_file_name, "wb"))
    , num_elts_in_buf_(0) {
    if (!output_file_) {
      fprintf(stderr, "Error: unable to open file. %s\n", output_file_name);
      exit(EXIT_FAILURE);
    }
  }

  ~ElementWriter() {
    if (is_buf_nonempty()) {
      fprintf(stderr, "Error: nonempty output buffer on completion. %i\n",
        num_elts_in_buf_);
      exit(EXIT_FAILURE);
    }
    fclose(output_file_);    
  }

  bool is_buf_full() const {
    return BUF_SIZE == num_elts_in_buf_;
  }

  bool is_buf_nonempty() const {
    return 0 != num_elts_in_buf_;
  }

  void add(int value) {

    // Output value has two bits.
    if (value < 0 || value > 3) {
      fprintf(stderr, "Invalid output value. %i\n", value);
      exit(EXIT_FAILURE);
    }

    // Internal consistency check.
    if (num_elts_in_buf_ < 0 || num_elts_in_buf_ >= BUF_SIZE) {
      fprintf(stderr, "Internal error: "
              "invalid number of elements in buffer. %i\n", num_elts_in_buf_);
      exit(EXIT_FAILURE);
    }

    // Add the two bits to the buffer.

    //buffer_ = (buffer_ << 2) | value;
    buffer_ = buffer_ | (value << 2*num_elts_in_buf_);;
    num_elts_in_buf_++;

    // Flush if full.

    if (is_buf_full())
      flush();
  }

  template<typename T> void add(std::string& value_string) const {


    std::stringstream value_stream(value_string);
    T value;
    value_stream >> value;
    const size_t num_written = fwrite(&value, sizeof(value), 1, output_file_);
    if (num_written != 1) {
      fprintf(stderr, "Error: file write error. %zu\n",
              num_written);
      exit(EXIT_FAILURE);
    }
  }

  void flush() {

    // Flush buffer to file.

    if (!num_elts_in_buf_)
      return;

    // NOTE: if the buffer is not full, then make low order
    // bits into high order bits to write out.

//    if (BUF_SIZE != num_elts_in_buf_)
//      buffer_ = buffer_ << (2 * (4 - num_elts_in_buf_));

    fputc(buffer_, output_file_);
    buffer_ = 0;
    num_elts_in_buf_ = 0;
  }

private:

  FILE* output_file_;
  int num_elts_in_buf_;
  char buffer_;
};

//-----------------------------------------------------------------------------

bool is_delim(int c) {
  return '\t' == c || ' ' == c;
}

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  // Parse arguments.

  if (argc != 4) {
    printf("preprocess: convert text SNP (tped) file "
           "to a binary SNP (CoMet input) file.\n");
    printf("Usage: preprocess <metric_type_prec> "
           "<snp_text_file> <snp_bin_file>\n");
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

  FILE* snptxtfile = fopen(argv[argnum], "r");
  if (!snptxtfile) {
    fprintf(stderr, "Error: unable to open file. %s\n", argv[argnum]);
    return 1;
  }
  argnum++;

  ElementWriter element_writer(argv[argnum]);
  argnum++;

  // Initializations.

  std::string line;

  const int num_frontmatter_fields = 4;

  int c = 0;
  int cprev = 0;
  int line_len = 0;
  int clo = 0;
  int chi = 0;
  int col = 0;

  //----------
  // LOOP to input chars from input file
  //----------

  while ((c = fgetc(snptxtfile)) != EOF) {

    if (line.size() < line_len+1)
      line.resize(line_len+1);

    // Store this character of line.
    line[line_len++] = c;

    // Wait for end of line.
    if (c != '\n')
      continue;

    // Reached end of a line, so process it.

    line_len--;

    //----------
    // PASS 1: loop over tokens to identify the allele labels.
    //----------

    // Only needed for ccc, duo.
    if (!is_czek) {

      clo = 0;
      chi = 0;

      for (int i=0, num_delims=0; i<line_len; ++i) {

        // Get character in line.
        c = line[i];
        if (is_delim(c)) {
          // Found delimiter, go to next token
          num_delims++;
          continue;
        }
        // Skip first four tokens ("frontmatter").
        if (num_delims < num_frontmatter_fields)
          continue;

        // Get token number
        col = num_delims - num_frontmatter_fields;

        // Perform check.
        if (col % 2 == 1 && is_ccc) {
          if ((c=='0' && cprev!='0') ||
              (c!='0' && cprev=='0')) {
            fprintf(stderr,
              "Error: token pair must be both zero or both nonzero\n");
            return 1;
          }
        }

//        // Handle Czekanowski case.
//
//        if (is_czek) {
//          // Skip to end of token.
//          while (i+1 < line_len && !is_delim(line[i+1])) {
//            ++i;
//          }
//
//          continue;
//        }

        // Record values of the tokens encountered.

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
      } // for i - PASS 1.

      //----------

      if (0 == chi)
        chi = clo;

      if (col % 2 == 0 && is_ccc) {
        fprintf(stderr, "Error: ccc requires an even number "
                "of non-frontmatter tokens. %i\n", col-1);
        return 1;
      }

    } // if (!is_czek)

    //----------
    // PASS 2: loop to output results.
    //----------

    for (int i=0, num_delims=0; i<line_len; ++i) {

      // Get character in line.
      c = line[i];
      if (is_delim(c)) {
        // Found delimiter, go to next token
        num_delims++;
        continue;
      }
      // Skip first four tokens ("frontmatter").
      if (num_delims < num_frontmatter_fields)
        continue;

      // Get token number.
      col = num_delims - num_frontmatter_fields;

      //----------
      // 1. Handle Czekanowski case.
      //----------

      if (is_czek) {
        // Find end of token.
        int len = 1;
        while (i+len < line_len && !is_delim(line[i+len])) {
          ++len;
        }

        // Write token to file in binary format.
        std::string token(line.substr(i, len));
        if (is_single)
          element_writer.add<float>(token);
        else
          element_writer.add<double>(token);

        // Complete the processing of this token.
        i += len - 1;

      } else {

      //----------
      // 2. Handle ccc, duo case.
      //----------

        // ccc case: only process when we have two input tokens.
        if (col % 2 == 0 && is_ccc) {
          cprev = c;
          continue;
        }

        const int c0 = cprev;
        const int c1 = c;

        enum {SN_UNDEFINED = -1};

        // Map the token or token pair to 2-bits element needed by CoMet.
        // NOTE the alphabetically first allele is mapped to the zero bit.

        // NOTE: order of lines matters below.
        int sn = is_duo ? (
              c1 == '0' ? 2 * (1) + 1 * (0)  //  0  ->  1 0
            : c1 == clo ? 2 * (0) + 1 * (0)  //  A  ->  0 0
            : c1 == chi ? 2 * (0) + 1 * (1)  //  B  ->  0 1
            : SN_UNDEFINED
          ) : (
              c0 == '0' && c1 == '0' ? 2 * (1) + 1 * (0)  //  00  ->  1 0
            : c0 == clo && c1 == clo ? 2 * (0) + 1 * (0)  //  AA  ->  0 0
            : c0 == clo && c1 == chi ? 2 * (0) + 1 * (1)  //  AB  ->  0 1
            : c0 == chi && c1 == clo ? 2 * (0) + 1 * (1)  //  BA  ->  0 1
            : c0 == chi && c1 == chi ? 2 * (1) + 1 * (1)  //  BB  ->  1 1
            : SN_UNDEFINED
          );

        if (SN_UNDEFINED == sn) {
          fprintf(stderr, "Error: error encountered mapping tokens.\n");
          return 1;
        }

        // Write 2 bits to file (buffered).
        element_writer.add(sn);

      } // if (is_czek)

      //----------
      // End handling of cases.
      //----------

    } // for i - PASS 2.

    //----------

    // Final flush of buffer

    if (element_writer.is_buf_nonempty())
      element_writer.flush();

    // Process next line

    line_len = 0;
    clo = 0;
    chi = 0;

  } // while c

  // Finish

  fclose(snptxtfile);

  return 0;
}

//-----------------------------------------------------------------------------

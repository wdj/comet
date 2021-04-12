
#include "stdio.h"

int main() {

  size_t v = 0;

  while (!feof (stdin)) {
    fscanf (stdin, "%zu", &v);
    fwrite(&v, 1, sizeof(v), stdout);
  }

}

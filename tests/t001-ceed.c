/// @file
/// Test creation and destruction of a CEED object
/// \test Test creation and destruction of a CEED object
#include <ceed.h>

int main(int argc, char **argv) {
  Ceed ceed;
  CeedInt type = -1;

  CeedInit(argv[1], &ceed);
  CeedGetPreferredMemType(ceed, (CeedMemType *)&type);

  if (type == -1)
    printf("Error getting preferred memory type. %d \n",type);

  CeedDestroy(&ceed);
  return 0;
}
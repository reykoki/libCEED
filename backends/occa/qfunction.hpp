// Copyright (c) 2019, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-734707.
// All Rights reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

#ifndef CEED_OCCA_QFUNCTION_HEADER
#define CEED_OCCA_QFUNCTION_HEADER

#include "ceed-object.hpp"
#include "qfunction-args.hpp"


namespace ceed {
  namespace occa {
    class QFunction : public CeedObject {
     public:
      // Ceed object information
      size_t ceedContextSize;
      void *ceedContext;

      // Owned resources
      std::string filename;
      std::string qFunctionName;
      ::occa::kernel qFunctionKernel;
      ::occa::memory qFunctionContext;
      QFunctionArgs args;

      QFunction(const std::string &source);

      ~QFunction();

      static QFunction* from(CeedQFunction qf);

      int buildKernel(const CeedInt Q);
      std::string getKernelSource(const std::string &kernelName,
                                  const CeedInt Q);

      void syncContext();

      int apply(CeedInt Q, CeedVector *U, CeedVector *V);

      //---[ Ceed Callbacks ]-----------
      static int registerQFunctionFunction(Ceed ceed, CeedQFunction qf,
                                           const char *fname, ceed::occa::ceedFunction f);

      static int ceedCreate(CeedQFunction qf);

      static int ceedApply(CeedQFunction qf,
                           CeedInt Q, CeedVector *U, CeedVector *V);

      static int ceedDestroy(CeedQFunction qf);
    };
  }
}

#endif

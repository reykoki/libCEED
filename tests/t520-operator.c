/// @file
/// Test creation, action, and destruction for composite mass matrix operator
/// \test Test creation, action, and destruction for composite mass matrix operator
#include <ceed.h>
#include <stdlib.h>
#include <math.h>
#include "t320-basis.h"
#include "t510-operator.h"

/* The mesh comprises of two rows of 3 quadralaterals followed by one row
     of 6 triangles:
   _ _ _
  |_|_|_|
  |_|_|_|
  |/|/|/|

*/

int main(int argc, char **argv) {
  Ceed ceed;
  CeedElemRestriction ErestrictxTet, ErestrictuTet,
                      ErestrictuiTet,
                      ErestrictxHex, ErestrictuHex,
                      ErestrictuiHex;
  CeedBasis bxTet, buTet,
            bxHex, buHex;
  CeedQFunction qf_setupTet, qf_massTet,
                qf_setupHex, qf_massHex;
  CeedOperator op_setupTet, op_massTet,
               op_setupHex, op_massHex,
               op_setup, op_mass;
  CeedVector qdataTet, qdataHex, X, U, V;
  const CeedScalar *hv;
  CeedInt nelemTet = 6, PTet = 6, QTet = 4,
          nelemHex = 6, PHex = 3, QHex = 4, dim = 2;
  CeedInt nx = 3, ny = 3,
          nxTet = 3, nyTet = 1, nxHex = 3;
  CeedInt row, col, offset;
  CeedInt ndofs = (nx*2+1)*(ny*2+1),
          nqptsTet = nelemTet*QTet,
          nqptsHex = nelemHex*QHex*QHex;
  CeedInt indxTet[nelemTet*PTet],
          indxHex[nelemHex*PHex*PHex];
  CeedScalar x[dim*ndofs];
  CeedScalar qref[dim*QTet], qweight[QTet];
  CeedScalar interp[PTet*QTet], grad[dim*PTet*QTet];

  CeedInit(argv[1], &ceed);

  // DoF Coordinates
  for (CeedInt i=0; i<ny*2+1; i++)
    for (CeedInt j=0; j<nx*2+1; j++) {
      x[i+j*(ny*2+1)+0*ndofs] = (CeedScalar) i / (2*ny);
      x[i+j*(ny*2+1)+1*ndofs] = (CeedScalar) j / (2*nx);
    }
  CeedVectorCreate(ceed, dim*ndofs, &X);
  CeedVectorSetArray(X, CEED_MEM_HOST, CEED_USE_POINTER, x);

  // Qdata Vectors
  CeedVectorCreate(ceed, nqptsTet, &qdataTet);
  CeedVectorCreate(ceed, nqptsHex, &qdataHex);

  // Set up Tet Elements
  for (CeedInt i=0; i<nelemTet/2; i++) {
    col = i % nxTet;
    row = i / nxTet;
    offset = col*2 + row*(nxTet*2+1)*2;

    indxTet[i*2*PTet +  0] =  2 + offset;
    indxTet[i*2*PTet +  1] =  9 + offset;
    indxTet[i*2*PTet +  2] = 16 + offset;
    indxTet[i*2*PTet +  3] =  1 + offset;
    indxTet[i*2*PTet +  4] =  8 + offset;
    indxTet[i*2*PTet +  5] =  0 + offset;

    indxTet[i*2*PTet +  6] = 14 + offset;
    indxTet[i*2*PTet +  7] =  7 + offset;
    indxTet[i*2*PTet +  8] =  0 + offset;
    indxTet[i*2*PTet +  9] = 15 + offset;
    indxTet[i*2*PTet + 10] =  8 + offset;
    indxTet[i*2*PTet + 11] = 16 + offset;
  }

  // -- Restrictions
  CeedElemRestrictionCreate(ceed, nelemTet, PTet, dim, ndofs, dim*ndofs,
                            CEED_MEM_HOST, CEED_USE_POINTER, indxTet,
                            &ErestrictxTet);

  CeedElemRestrictionCreate(ceed, nelemTet, PTet, 1, 1, ndofs,
                            CEED_MEM_HOST, CEED_USE_POINTER, indxTet,
                            &ErestrictuTet);
  CeedInt stridesuTet[3] = {1, QTet, QTet};
  CeedElemRestrictionCreateStrided(ceed,  nelemTet, QTet, 1, nqptsTet,
                                   stridesuTet, &ErestrictuiTet);

  // -- Bases
  buildmats(qref, qweight, interp, grad);
  CeedBasisCreateH1(ceed, CEED_TRIANGLE, dim, PTet, QTet, interp, grad, qref,
                    qweight, &bxTet);

  buildmats(qref, qweight, interp, grad);
  CeedBasisCreateH1(ceed, CEED_TRIANGLE, 1, PTet, QTet, interp, grad, qref,
                    qweight, &buTet);

  // -- QFunctions
  CeedQFunctionCreateInterior(ceed, 1, setup, setup_loc, &qf_setupTet);
  CeedQFunctionAddInput(qf_setupTet, "_weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddInput(qf_setupTet, "dx", dim*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddOutput(qf_setupTet, "rho", 1, CEED_EVAL_NONE);

  CeedQFunctionCreateInterior(ceed, 1, mass, mass_loc, &qf_massTet);
  CeedQFunctionAddInput(qf_massTet, "rho", 1, CEED_EVAL_NONE);
  CeedQFunctionAddInput(qf_massTet, "u", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(qf_massTet, "v", 1, CEED_EVAL_INTERP);

  // -- Operators
  // ---- Setup Tet
  CeedOperatorCreate(ceed, qf_setupTet, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE, &op_setupTet);
  CeedOperatorSetField(op_setupTet, "_weight", CEED_ELEMRESTRICTION_NONE, bxTet,
                       CEED_VECTOR_NONE);
  CeedOperatorSetField(op_setupTet, "dx", ErestrictxTet, bxTet,
                       CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_setupTet, "rho", ErestrictuiTet,
                       CEED_BASIS_COLLOCATED, qdataTet);
  // ---- Mass Tet
  CeedOperatorCreate(ceed, qf_massTet, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                     &op_massTet);
  CeedOperatorSetField(op_massTet, "rho", ErestrictuiTet, CEED_BASIS_COLLOCATED,
                       qdataTet);
  CeedOperatorSetField(op_massTet, "u", ErestrictuTet, buTet,
                       CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_massTet, "v", ErestrictuTet, buTet,
                       CEED_VECTOR_ACTIVE);

  // Set up Hex Elements
  for (CeedInt i=0; i<nelemHex; i++) {
    col = i % nxHex;
    row = i / nxHex;
    offset = (nxTet*2+1)*(nyTet*2)*(1+row) + col*2;
    for (CeedInt j=0; j<PHex; j++)
      for (CeedInt k=0; k<PHex; k++)
        indxHex[PHex*(PHex*i+k)+j] = offset + k*(nxHex*2+1) + j;
  }

  // -- Restrictions
  CeedElemRestrictionCreate(ceed, nelemHex, PHex*PHex, dim, ndofs, dim*ndofs,
                            CEED_MEM_HOST, CEED_USE_POINTER, indxHex,
                            &ErestrictxHex);

  CeedElemRestrictionCreate(ceed, nelemHex, PHex*PHex, 1, 1, ndofs,
                            CEED_MEM_HOST, CEED_USE_POINTER, indxHex,
                            &ErestrictuHex);
  CeedInt stridesuHex[3] = {1, QHex*QHex, QHex*QHex};
  CeedElemRestrictionCreateStrided(ceed, nelemHex, QHex*QHex, 1, nqptsHex,
                                   stridesuHex, &ErestrictuiHex);

  // -- Bases
  CeedBasisCreateTensorH1Lagrange(ceed, dim, dim, PHex, QHex, CEED_GAUSS,
                                  &bxHex);
  CeedBasisCreateTensorH1Lagrange(ceed, dim, 1, PHex, QHex, CEED_GAUSS, &buHex);

  // -- QFunctions
  CeedQFunctionCreateInterior(ceed, 1, setup, setup_loc, &qf_setupHex);
  CeedQFunctionAddInput(qf_setupHex, "_weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddInput(qf_setupHex, "dx", dim*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddOutput(qf_setupHex, "rho", 1, CEED_EVAL_NONE);

  CeedQFunctionCreateInterior(ceed, 1, mass, mass_loc, &qf_massHex);
  CeedQFunctionAddInput(qf_massHex, "rho", 1, CEED_EVAL_NONE);
  CeedQFunctionAddInput(qf_massHex, "u", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(qf_massHex, "v", 1, CEED_EVAL_INTERP);

  // -- Operators
  CeedOperatorCreate(ceed, qf_setupHex, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE, &op_setupHex);
  CeedOperatorSetField(op_setupHex, "_weight", CEED_ELEMRESTRICTION_NONE, bxHex,
                       CEED_VECTOR_NONE);
  CeedOperatorSetField(op_setupHex, "dx", ErestrictxHex, bxHex,
                       CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_setupHex, "rho", ErestrictuiHex,
                       CEED_BASIS_COLLOCATED, qdataHex);

  CeedOperatorCreate(ceed, qf_massHex, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                     &op_massHex);
  CeedOperatorSetField(op_massHex, "rho", ErestrictuiHex, CEED_BASIS_COLLOCATED,
                       qdataHex);
  CeedOperatorSetField(op_massHex, "u", ErestrictuHex, buHex,
                       CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_massHex, "v", ErestrictuHex, buHex,
                       CEED_VECTOR_ACTIVE);

  // Set up Composite Operators
  // -- Create
  CeedCompositeOperatorCreate(ceed, &op_setup);
  // -- Add SubOperators
  CeedCompositeOperatorAddSub(op_setup, op_setupTet);
  CeedCompositeOperatorAddSub(op_setup, op_setupHex);

  // -- Create
  CeedCompositeOperatorCreate(ceed, &op_mass);
  // -- Add SubOperators
  CeedCompositeOperatorAddSub(op_mass, op_massTet);
  CeedCompositeOperatorAddSub(op_mass, op_massHex);

  // Apply Setup Operator
  CeedOperatorApply(op_setup, X, CEED_VECTOR_NONE, CEED_REQUEST_IMMEDIATE);

  // Apply Mass Operator
  CeedVectorCreate(ceed, ndofs, &U);
  CeedVectorSetValue(U, 0.0);
  CeedVectorCreate(ceed, ndofs, &V);

  CeedOperatorApply(op_mass, U, V, CEED_REQUEST_IMMEDIATE);

  // Check output
  CeedVectorGetArrayRead(V, CEED_MEM_HOST, &hv);
  for (CeedInt i=0; i<ndofs; i++)
    if (fabs(hv[i]) > 1e-14) printf("[%d] v %g != 0.0\n",i, hv[i]);
  CeedVectorRestoreArrayRead(V, &hv);

  // Cleanup
  CeedQFunctionDestroy(&qf_setupTet);
  CeedQFunctionDestroy(&qf_massTet);
  CeedOperatorDestroy(&op_setupTet);
  CeedOperatorDestroy(&op_massTet);
  CeedQFunctionDestroy(&qf_setupHex);
  CeedQFunctionDestroy(&qf_massHex);
  CeedOperatorDestroy(&op_setupHex);
  CeedOperatorDestroy(&op_massHex);
  CeedOperatorDestroy(&op_setup);
  CeedOperatorDestroy(&op_mass);
  CeedElemRestrictionDestroy(&ErestrictuTet);
  CeedElemRestrictionDestroy(&ErestrictxTet);
  CeedElemRestrictionDestroy(&ErestrictuiTet);
  CeedElemRestrictionDestroy(&ErestrictuHex);
  CeedElemRestrictionDestroy(&ErestrictxHex);
  CeedElemRestrictionDestroy(&ErestrictuiHex);
  CeedBasisDestroy(&buTet);
  CeedBasisDestroy(&bxTet);
  CeedBasisDestroy(&buHex);
  CeedBasisDestroy(&bxHex);
  CeedVectorDestroy(&X);
  CeedVectorDestroy(&U);
  CeedVectorDestroy(&V);
  CeedVectorDestroy(&qdataTet);
  CeedVectorDestroy(&qdataHex);
  CeedDestroy(&ceed);
  return 0;
}

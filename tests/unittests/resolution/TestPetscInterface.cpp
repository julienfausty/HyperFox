#include <catch2/catch.hpp>
#include <string>
#include <map>
#include <tuple>
#include "PetscInterface.h"
#include "petscksp.h"
#include "PetscOpts.h"

using namespace hfox;

TEST_CASE("Testing the PetscInterface", "[unit][resolution][PetscInterface]"){

  SECTION("Construction"){
    CHECK_NOTHROW(PetscInterface());
    PetscOpts pOpts;
    pOpts.rtol = 1e-10;
    pOpts.solverType = KSPCG;
    pOpts.verbose=0;
    CHECK_NOTHROW(PetscInterface(pOpts));
  };

  SECTION("Test setVal"){
    double val = 1.0;
    double v;
    int i = 1, j = 2, n = 3;
    PetscInterface pIFace;
    pIFace.initialize(); pIFace.configure(); pIFace.allocate(n);
    CHECK_NOTHROW(pIFace.setValMatrix(i, j, val));
    CHECK_NOTHROW(pIFace.setValRHS(i, val));
    CHECK_NOTHROW(pIFace.assemble());
    const Mat * M = pIFace.getMatrix();
    const Vec * rhs = pIFace.getRHS();
    MatGetValues(*M, 1, &i, 1, &j, &v);
    CHECK(v == val);
    VecGetValues(*rhs, 1, &i, &v);
    CHECK(v == val);
  };


  SECTION("Test setVals"){
    int n = 4;
    std::vector<double> vals = 
    {1.0, 0.0, 0.0, 0.0, 
      0.0, 2.0, 0.0, 0.0,
      0.0, 0.0, 3.0, 0.0, 
      0.0, 0.0, 0.0, 4.0};
    double v;
    std::vector<int> is = {0, 1, 2, 3};
    std::vector<int> js = {0, 1, 2, 3};
    PetscInterface pIFace;
    pIFace.initialize(); pIFace.configure(); pIFace.allocate(n);
    CHECK_NOTHROW(pIFace.setValsMatrix(is, js, vals.data()));
    CHECK_NOTHROW(pIFace.setValsRHS(is, vals.data()));
    CHECK_NOTHROW(pIFace.assemble());
    const Mat * M = pIFace.getMatrix();
    const Vec * rhs = pIFace.getRHS();
    for(int i = 0; i < is.size(); i++){
      for(int j = 0; j < js.size(); j++){
        MatGetValues(*M, 1, &(is[i]), 1, &(js[j]), &v);
        CHECK(v == vals[i*is.size() + j]);
      }
      VecGetValues(*rhs, 1, &(is[i]), &v);
      CHECK(v == vals[i]);
    }
  };


  SECTION("Test addVal"){
    double val = 1.0;
    double v;
    int i = 1, j = 2, n = 3;
    PetscInterface pIFace;
    pIFace.initialize(); pIFace.configure(); pIFace.allocate(n);
    CHECK_NOTHROW(pIFace.addValMatrix(i, j, val));
    CHECK_NOTHROW(pIFace.addValRHS(i, val));
    CHECK_NOTHROW(pIFace.assemble());
    const Mat * M = pIFace.getMatrix();
    const Vec * rhs = pIFace.getRHS();
    MatGetValues(*M, 1, &i, 1, &j, &v);
    CHECK(v == val);
    VecGetValues(*rhs, 1, &i, &v);
    CHECK(v == val);
  };


  SECTION("Test addVals"){
    int n = 4;
    std::vector<double> vals = 
    {1.0, 0.0, 0.0, 0.0, 
      0.0, 2.0, 0.0, 0.0,
      0.0, 0.0, 3.0, 0.0, 
      0.0, 0.0, 0.0, 4.0};
    std::vector<double> skew = 
    {0.0, 2.0, 0.0, 0.0, 
      -2.0, 0.0, 1.0, 0.0,
      0.0, -1.0, 0.0, -3.0, 
      0.0, 0.0, 3.0, 0.0};
    double v;
    std::vector<int> is = {0, 1, 2, 3};
    std::vector<int> js = {0, 1, 2, 3};
    PetscInterface pIFace;
    pIFace.initialize(); pIFace.configure(); pIFace.allocate(n);
    CHECK_NOTHROW(pIFace.addValsMatrix(is, js, vals.data()));
    CHECK_NOTHROW(pIFace.addValsMatrix(is, js, skew.data()));
    CHECK_NOTHROW(pIFace.addValsRHS(is, vals.data()));
    CHECK_NOTHROW(pIFace.assemble());
    const Mat * M = pIFace.getMatrix();
    const Vec * rhs = pIFace.getRHS();
    for(int i = 0; i < is.size(); i++){
      for(int j = 0; j < js.size(); j++){
        MatGetValues(*M, 1, &(is[i]), 1, &(js[j]), &v);
        CHECK(v == vals[i*is.size() + j] + skew[i*is.size() + j]);
      }
      VecGetValues(*rhs, 1, &(is[i]), &v);
      CHECK(v == vals[i]);
    }
  };

};

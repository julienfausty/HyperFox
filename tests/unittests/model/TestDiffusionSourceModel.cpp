#include <catch2/catch.hpp>
#include <string>
#include "DiffusionSource.h"
#include "Euler.h"

using namespace hfox;

double testSrcDiffSrc(const std::vector<double> & v){
  return v[0];
};

TEST_CASE("Testing the DiffusionSource model", "[unit][model][DiffusionSource]"){
  
  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      SECTION("Testing setup (dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(DiffusionSource(&refEl));
        DiffusionSource mod(&refEl);
        CHECK_THROWS(mod.compute());
        std::map<std::string, std::vector<double> > fm;
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod.setFieldMap(&fm));
        CHECK_THROWS(mod.compute());
        CHECK_THROWS(mod.setSourceFunction(testSrcDiffSrc));
        CHECK_NOTHROW(mod.allocate(1));
        CHECK_NOTHROW(mod.setSourceFunction(testSrcDiffSrc));
      };

      DiffusionSource mod(&refEl);
      std::map<std::string, std::vector<double> > fm;
      std::vector<double> sol(refEl.getNumNodes(), 0.0);
      std::iota(sol.begin(), sol.end(), 0.0);
      fm["Solution"] = sol;
      mod.setFieldMap(&fm);
      mod.allocate(1);
      mod.setSourceFunction(testSrcDiffSrc);
      Diffusion diffOp(&refEl);
      Source srcOp(&refEl);
      diffOp.allocate(1);
      srcOp.allocate(1);
      srcOp.setSourceFunction(testSrcDiffSrc);
      std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
      std::vector<double> dV(refEl.getNumIPs(), 0.0);
      std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
      diffOp.assemble(dV, invJacs);
      srcOp.calcSource(*(refEl.getNodes()));
      srcOp.assemble(dV, invJacs);
      SECTION("Poisson equation in reference element"){
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod.compute());
        EMatrix testMat = *(diffOp.getMatrix());
        testMat -= *(mod.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = (*(mod.getLocalRHS()));
        testVec -= ((EVector) *(srcOp.getMatrix()));
        CHECK(testVec.dot(testVec) < tol);
      };
      DiffusionSource mod2(&refEl);
      std::vector<double> diff(refEl.getNumNodes(), 3.0);
      fm["DiffusionTensor"] = diff;
      Euler ts(&refEl);
      double timeStep = 1e-2;
      ts.setTimeStep(timeStep);
      Mass mass(&refEl);
      mass.allocate(1);
      mass.assemble(dV, invJacs);
      mod2.setTimeScheme(&ts);
      mod2.setFieldMap(&fm);
      mod2.allocate(1);
      mod2.setSourceFunction(testSrcDiffSrc);
      SECTION("Diffusion source equation in reference element"){
        CHECK_NOTHROW(mod2.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod2.compute());
        EMatrix testMat = (*(diffOp.getMatrix()));
        testMat *= 3*timeStep;
        testMat += (*(mass.getMatrix()));
        testMat -= *(mod2.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec(testMat.rows());
        testVec = timeStep*((EVector) srcOp.getMatrix()->block(0, 0, testMat.rows(), 1));
        testVec += (*(mass.getMatrix()))*EMap<EVector>(sol.data(), sol.size());
        testVec -= (*(mod2.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };
    }
  }

};

#include <catch2/catch.hpp>
#include "BurgersModel.h"
#include "Euler.h"

using namespace hfox;

TEST_CASE("Testing the BurgersModel", "[unit][model][nonlinear][BurgersModel]"){
  
  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      SECTION("Testing setup (dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(BurgersModel(&refEl));
        BurgersModel mod(&refEl);
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.allocate(i+1));
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        std::map<std::string, std::vector<double> > fm;
        CHECK_THROWS(mod.setFieldMap(&fm));
        fm["BufferSolution"] = std::vector<double>(refEl.getNumNodes()*(i+1), 1.0);
        CHECK_NOTHROW(mod.setFieldMap(&fm));
      };

      BurgersModel mod(&refEl);
      std::map<std::string, std::vector<double> > fm;
      std::vector<double> sol(refEl.getNumNodes() * (i+1), 0.0);
      std::vector<double> buffSol(refEl.getNumNodes() * (i+1), 2.0);
      std::vector<double> diff(refEl.getNumNodes()*(i+1)*(i+1), 0);
      for(int k = 0; k < refEl.getNumNodes(); k++){
        EMap<EMatrix>(diff.data() + (i+1)*(i+1)*k, (i+1), (i+1)) = EMatrix::Identity((i+1), (i+1))*3.0;
      }
      std::iota(sol.begin(), sol.end(), 0.0);
      fm["Solution"] = sol;
      fm["BufferSolution"] = buffSol;
      fm["DiffusionTensor"] = diff;
      mod.setFieldMap(&fm);
      mod.allocate(i+1);
      std::vector<EVector> velVectors(refEl.getNumNodes(), EVector::Constant(i+1, 2.0));
      NabUU nabuu(&refEl);
      std::vector<EMatrix> diffTensors(refEl.getNumNodes(), EMatrix::Identity(i+1, i+1)*3.0);
      Diffusion diffOp(&refEl);
      nabuu.allocate(i+1);
      nabuu.setSolution(velVectors);
      diffOp.allocate(i+1);
      diffOp.setDiffTensor(diffTensors);
      std::vector<EMatrix> jacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
      std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
      std::vector<double> dV(refEl.getNumIPs(), 0.0);
      std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
      nabuu.assemble(dV, invJacs);
      diffOp.assemble(dV, invJacs);
      SECTION("No time operator in reference element(dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod.compute());
        EMatrix testMat = *(nabuu.getMatrix()) + *(diffOp.getMatrix());
        testMat -= *(mod.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = (*(mod.getLocalRHS()));
        testVec -= *(nabuu.getMatrix()) * EMap<EVector>(buffSol.data(), buffSol.size());
        CHECK(testVec.dot(testVec) < tol);
      };

      BurgersModel mod2(&refEl);
      Euler ts(&refEl);
      double timeStep = 1e-2;
      ts.setTimeStep(timeStep);
      Mass mass(&refEl);
      mass.allocate(i+1);
      mass.assemble(dV, invJacs);
      mod2.setTimeScheme(&ts);
      mod2.allocate(i+1);
      mod2.setFieldMap(&fm);
      SECTION("Transport equation in reference element(dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(mod2.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod2.compute());
        EMatrix testMat = *(nabuu.getMatrix()) + *(diffOp.getMatrix());
        testMat *= timeStep;
        testMat += (*(mass.getMatrix()));
        testMat -= *(mod2.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = *(nabuu.getMatrix()) * EMap<EVector>(buffSol.data(), buffSol.size()) * timeStep;
        testVec = (*(mass.getMatrix()))*EMap<EVector>(sol.data(), sol.size());
        testVec -= (*(mod2.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };
    }
  }
};
#include <catch2/catch.hpp>
#include <string>
#include "Transport.h"
#include "Convection.h"
#include "Mass.h"
#include "Euler.h"

using namespace hfox;

TEST_CASE("Testing Transport model", "[unit][model][Transport]"){

  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 0; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      SECTION("Testing setup (dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(Transport(&refEl));
        Transport mod(&refEl);
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.allocate(1));
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        std::map<std::string, std::vector<double> > fm;
        CHECK_THROWS(mod.setFieldMap(&fm));
        fm["Velocity"] = std::vector<double>(refEl.getNumNodes()*(i+1), 1.0);
        CHECK_NOTHROW(mod.setFieldMap(&fm));
      };
      
      Transport mod(&refEl);
      std::map<std::string, std::vector<double> > fm;
      std::vector<double> sol(refEl.getNumNodes(), 0.0);
      std::iota(sol.begin(), sol.end(), 0.0);
      fm["Solution"] = sol;
      double v = 2.0;
      std::vector<double> vel(refEl.getNumNodes()*(i+1), 0.0);
      std::vector<EVector> velVect(refEl.getNumNodes(), EVector(i+1));
      for(int k = 0; k < refEl.getNumNodes(); k++){
        vel[k*(i+1)] = v;
        velVect[k][0] = v;
      }
      fm["Velocity"] = vel;
      mod.setFieldMap(&fm);
      mod.allocate(1);
      Convection convOp(&refEl);
      convOp.allocate(1);
      std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
      std::vector<double> dV(refEl.getNumIPs(), 0.0);
      std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
      convOp.setVelocity(velVect);
      convOp.assemble(dV, invJacs);
      SECTION("Grad ortho equation in reference element(dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod.compute());
        EMatrix testMat = *(convOp.getMatrix());
        testMat -= *(mod.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = (*(mod.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };
      Transport mod2(&refEl);
      Euler ts(&refEl);
      double timeStep = 1e-2;
      ts.setTimeStep(timeStep);
      Mass mass(&refEl);
      mass.allocate(1);
      mass.assemble(dV, invJacs);
      mod2.setTimeScheme(&ts);
      mod2.allocate(1);
      mod2.setFieldMap(&fm);
      SECTION("Transport equation in reference element(dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(mod2.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod2.compute());
        EMatrix testMat = timeStep*(*(convOp.getMatrix()));
        testMat += (*(mass.getMatrix()));
        testMat -= *(mod2.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec(testMat.rows());
        testVec = (*(mass.getMatrix()))*EMap<EVector>(sol.data(), sol.size());
        testVec -= (*(mod2.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };
    }
  }
};

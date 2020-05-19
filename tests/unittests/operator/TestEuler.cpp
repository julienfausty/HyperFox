#include <catch2/catch.hpp>
#include "Euler.h"
#include "Mass.h"
#include "Convection.h"

using namespace hfox;

TEST_CASE("Testing the Euler time scheme", "[unit][operator][TimeScheme][Euler]"){

  SECTION("Testing necessary"){
    ReferenceElement refEl(2, 1, "simplex");
    CHECK_NOTHROW(Euler(&refEl, false));
    CHECK_NOTHROW(Euler(&refEl, true));
    CHECK_NOTHROW(Euler(&refEl));
    Euler backEuler(&refEl);
    EMatrix * dummyMat;
    EVector * dummyVec;
    CHECK_THROWS(backEuler.apply(dummyMat, dummyVec));
    std::map<std::string, std::vector<double> > fm;
    CHECK_THROWS(backEuler.setFieldMap(&fm));
    std::vector<double> sol(refEl.getNumNodes(), 1.0);
    fm["Solution"] = sol;
    CHECK_NOTHROW(backEuler.setFieldMap(&fm));
  };

  double tol = 1e-12;

  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 0; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      SECTION("Testing in reference element (dim= " + std::to_string(i+1) + ", order= " + std::to_string(j+1) + ")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        Mass mass(&refEl);
        Convection conv(&refEl);
        Euler backEuler(&refEl);
        Euler forEuler(&refEl, true);
        mass.allocate(1);
        conv.allocate(1);
        backEuler.allocate(1);
        forEuler.allocate(1);
        EVector vel = EVector::Constant(i+1, 2.0);
        std::vector<EVector> velocity(refEl.getNumNodes(), vel);
        conv.setVelocity(velocity);
        EMatrix buffMat(refEl.getNumNodes(), refEl.getNumNodes());
        EMatrix anaMat = buffMat;
        EVector buffVec(refEl.getNumNodes());
        EVector anaVec = buffMat;
        std::vector<double> dV(refEl.getNumIPs(), 0.0);
        std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
        std::vector<EMatrix> jacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
        mass.assemble(dV, jacs);
        conv.assemble(dV, jacs);
        backEuler.assemble(dV, jacs);
        forEuler.assemble(dV, jacs);
        std::map<std::string, std::vector<double> > fm;
        std::vector<double> sol(refEl.getNumNodes(), 0.0);
        std::iota(sol.begin(), sol.end(), 0.0);
        fm["Solution"] = sol;
        backEuler.setFieldMap(&fm);
        forEuler.setFieldMap(&fm);
        double tStep = 1e-1;
        backEuler.setTimeStep(tStep);
        forEuler.setTimeStep(tStep);
        buffMat = *(conv.getMatrix());
        buffVec = EVector::Constant(refEl.getNumNodes(), 1.0);
        backEuler.apply(&buffMat, &buffVec);
        anaMat = (*(mass.getMatrix())) + tStep*(*(conv.getMatrix()));
        anaVec = (*(mass.getMatrix()))*EMap<const EVector>(sol.data(), sol.size());
        anaVec += tStep*EVector::Constant(refEl.getNumNodes(), 1.0);
        EMatrix testMat = anaMat-buffMat;
        EVector testVec = anaVec-buffVec;
        CHECK((testMat.transpose()*testMat).sum() < tol);
        CHECK(testVec.dot(testVec) < tol);
        buffMat = *(conv.getMatrix());
        buffVec = EVector::Constant(refEl.getNumNodes(), 1.0);
        forEuler.apply(&buffMat, &buffVec);
        anaMat = (*(mass.getMatrix()));
        anaVec = ((*(mass.getMatrix())) - tStep*(*(conv.getMatrix())))*EMap<const EVector>(sol.data(), sol.size());        
        anaVec += tStep*EVector::Constant(refEl.getNumNodes(), 1.0);
        testMat = anaMat-buffMat;
        testVec = anaVec-buffVec;
        CHECK((testMat.transpose()*testMat).sum() < tol);
        CHECK(testVec.dot(testVec) < tol);
      };
    }
  }
}

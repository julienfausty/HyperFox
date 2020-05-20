#include <catch2/catch.hpp>
#include "RungeKutta.h"
#include "Field.h"
#include "Convection.h"

using namespace hfox;

TEST_CASE("Testing RungeKutta class", "[unit][operator][TimeScheme][RungeKutta]"){

  SECTION("Testing necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    CHECK_NOTHROW(RungeKutta(&refEl));
    RKType t = FEuler;
    CHECK_NOTHROW(RungeKutta(&refEl, t));
    RungeKutta cn(&refEl);
    EMatrix * dummyMat;
    EVector * dummyVec;
    CHECK_THROWS(cn.apply(dummyMat, dummyVec));
    std::map<std::string, std::vector<double> > fm;
    CHECK_THROWS(cn.setFieldMap(&fm));
    std::vector<double> sol(refEl.getNumNodes(), 1.0);
    fm["Solution"] = sol;
    CHECK_THROWS(cn.setFieldMap(&fm));
    fm["OldSolution"] = sol;
    CHECK_THROWS(cn.setFieldMap(&fm));
    fm["RKStage_0"] = sol;
    CHECK_THROWS(cn.setFieldMap(&fm));
    fm["RKStage_1"] = sol;
    CHECK_NOTHROW(cn.setFieldMap(&fm));
  };

  SECTION("Testing Butcher table characteristics"){
    ReferenceElement refEl(1, 1, "simplex");
    RungeKutta ts(&refEl);
    CHECK(ts.getNumStages() == 2);
    CHECK(ts.getStage() == 0);
    //check upper triangle vals throw
    CHECK_THROWS(ts.setButcherTable(EMatrix::Constant(2, 2, 1)));
    //check non square matrix
    CHECK_THROWS(ts.setButcherTable(EMatrix::Identity(3, 3).block(0, 0, 3, 2)));
    CHECK_NOTHROW(ts.setButcherTable(EMatrix::Identity(3,3)));
    CHECK_NOTHROW(ts.setButcherTable(RK4));
    CHECK(ts.getNumStages() == 4);
  };

  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 0; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      SECTION("Testing Crank Nicolson (dim="+std::to_string(i+1)+", order="+std::to_string(j+1)+")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        RungeKutta ts(&refEl);
        Mass mass(&refEl);
        Convection conv(&refEl);
        mass.allocate(1);
        conv.allocate(1);
        ts.allocate(1);
        EVector vel = EVector::Constant(i+1, 2.0);
        std::vector<EVector> velocity(refEl.getNumNodes(), vel);
        conv.setVelocity(velocity);
        std::vector<double> dV(refEl.getNumIPs(), 0.0);
        std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
        std::vector<EMatrix> jacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
        mass.assemble(dV, jacs);
        conv.assemble(dV, jacs);
        ts.assemble(dV, jacs);
        std::map<std::string, std::vector<double> > fm;
        std::vector<double> sol(refEl.getNumNodes(), 0.0);
        std::iota(sol.begin(), sol.end(), 0.0);
        fm["Solution"] = sol;
        fm["OldSolution"] = sol;
        for(int k = 0; k < ts.getNumStages(); k++){
          fm["RKStage_"+std::to_string(k)] = sol;
        }
        ts.setFieldMap(&fm);
        sol.resize(0);
        sol.resize(refEl.getNumNodes(), 1.0);
        std::map<std::string, Field* > fieldMap;
        Field solF; *(solF.getValues()) = sol;
        fieldMap["Solution"] = &solF;
        sol.resize(0);
        sol.resize(refEl.getNumNodes(), 2.0);
        Field oldSolF; *(oldSolF.getValues()) = sol;
        fieldMap["OldSolution"] = &oldSolF;
        std::vector<Field*> rkStages(ts.getNumStages());
        for(int k = 0; k < ts.getNumStages(); k++){
          sol.resize(0);
          sol.resize(refEl.getNumNodes(), 3.0+k);
          rkStages[k] = new Field();
          *(rkStages[k]->getValues()) = sol;
          fieldMap["RKStage_"+std::to_string(k)] = rkStages[k];
        }
        double tStep = 1e-1;
        ts.setTimeStep(tStep);
        EMatrix buffMat = EMatrix::Zero(refEl.getNumNodes(), refEl.getNumNodes());
        EMatrix anaMat = buffMat;
        EVector buffVec = EVector::Zero(refEl.getNumNodes());
        EVector anaVec = buffVec;
        //0th step
        buffMat = *(conv.getMatrix());
        buffVec = EVector::Constant(refEl.getNumNodes(), 1.0);
        ts.apply(&buffMat, &buffVec);
        anaMat = (*(mass.getMatrix()));
        anaVec = ((*(mass.getMatrix())) - tStep*(*(conv.getMatrix())))*EMap<const EVector>(fm["OldSolution"].data(), refEl.getNumNodes());        
        anaVec += tStep*EVector::Constant(refEl.getNumNodes(), 1.0);
        EMatrix testMat = anaMat-buffMat;
        EVector testVec = anaVec-buffVec;
        CHECK((testMat.transpose()*testMat).sum() < tol);
        CHECK(testVec.dot(testVec) < tol);
        CHECK(ts.getStage() == 0);
        ts.computeStage(&fieldMap);
        CHECK(ts.getStage() == 1);
        double RKStageVal = (1.0-2.0)/tStep;
        for(int k = 0; k < refEl.getNumNodes(); k++){
          CHECK(fieldMap["RKStage_0"]->getValues()->at(k) == RKStageVal);
        }
        //1st step
        buffMat = *(conv.getMatrix());
        buffVec = EVector::Constant(refEl.getNumNodes(), 1.0);
        ts.apply(&buffMat, &buffVec);
        anaMat = (*(mass.getMatrix())) + (tStep/2.0)*(*(conv.getMatrix()));
        anaVec = EMap<const EVector>(fm["OldSolution"].data(), refEl.getNumNodes()) + (tStep/2.0)*(EMap<const EVector>(fm["RKStage_0"].data(), refEl.getNumNodes()));
        anaVec = (*(mass.getMatrix()))*EMap<const EVector>(fm["OldSolution"].data(), refEl.getNumNodes()) - tStep*(*(conv.getMatrix()))*anaVec;
        anaVec += tStep*EVector::Constant(refEl.getNumNodes(), 1.0);
        testMat = anaMat-buffMat;
        testVec = anaVec-buffVec;
        CHECK((testMat.transpose()*testMat).sum() < tol);
        CHECK(testVec.dot(testVec) < tol);
        CHECK(ts.getStage() == 1);
        ts.computeStage(&fieldMap);
        CHECK(ts.getStage() == 2);
        for(int k = 0; k < refEl.getNumNodes(); k++){
          CHECK(fieldMap["RKStage_1"]->getValues()->at(k) == RKStageVal);
        }
        //compute solution
        ts.computeSolution(&fieldMap);
        CHECK(ts.getStage() == 0);
        for(int k = 0; k < refEl.getNumNodes(); k++){
          CHECK(fieldMap["Solution"]->getValues()->at(k) == (fieldMap["OldSolution"]->getValues()->at(k) + tStep*RKStageVal));
        }
        for(int k = 0; k < ts.getNumStages(); k++){
          delete rkStages[k];
        }
      };
    }
  }
};

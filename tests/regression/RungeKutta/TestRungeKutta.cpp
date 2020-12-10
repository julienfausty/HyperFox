#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>
#include <fstream>
#include "Mesh.h"
#include "Mass.h"
#include "RungeKutta.h"
#include "TestUtils.h"

using namespace hfox;

namespace rungekutta{

typedef std::tuple<RKType, double, double> SimRes; 

void runRungeKutta(SimRes * res){
  //init reference element
  ReferenceElement refEl(1, 1, "simplex");
  //init ref element related quantities
  std::vector<double> dV = *(refEl.getIPWeights());
  std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(1, 1));
  //init time scheme
  RungeKutta ts(&refEl, std::get<0>(*res));
  int nStages = ts.getNumStages();
  //init fields
  Field u;
  *(u.getNumEntities()) = 1;
  *(u.getNumObjPerEnt()) = refEl.getNumNodes();
  *(u.getNumValsPerObj()) = 1;
  u.getValues()->resize(refEl.getNumNodes(), 0.0);
  Field oldU = u;
  std::vector<Field> rkStages(nStages, u);
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &u;
  fieldMap["OldSolution"] = &oldU;
  for(int i = 0; i < nStages; i++){
    fieldMap["RKStage_"+std::to_string(i)] = &(rkStages[i]);
  }
  //init fields
  std::fill(u.getValues()->begin(), u.getValues()->end(), 1.0);
  //init mass
  Mass m(&refEl);
  //allocate ts and m
  ts.allocate(1);
  m.allocate(1);
  //assemble local matrixes
  ts.assemble(dV, invJacs);
  m.assemble(dV, invJacs);
  double t = 0.0;
  double tend = 1.0;
  double tstep = std::get<1>(*res);
  ts.setTimeStep(tstep);
  int nIters = tend/tstep;
  double err2 = 0.0;
  for(int iter = 0; iter < nIters; iter++){
    t += tstep;
    //field updating
    std::copy(u.getValues()->begin(), u.getValues()->end(), oldU.getValues()->begin());
    for(int i = 0; i < nStages; i++){
      std::fill(rkStages[i].getValues()->begin(), rkStages[i].getValues()->end(), 0.0);
    }
    for(int i = 0; i < nStages; i++){
      //local fields
      std::map<std::string, std::vector<double> > localFieldMap;
      localFieldMap["Solution"] = *(u.getValues());
      localFieldMap["OldSolution"] = *(oldU.getValues());
      for(int i = 0; i < nStages; i++){
        localFieldMap["RKStage_"+std::to_string(i)] = *(rkStages[i].getValues());
      }
      //setup
      EMatrix stiff = -*(m.getMatrix());
      EVector rhs = EVector::Constant(refEl.getNumNodes(), 0.0);
      ts.setFieldMap(&localFieldMap);
      ts.apply(&stiff, &rhs);
      //solve
      EMap<EVector>(u.getValues()->data(), u.getValues()->size()) = stiff.colPivHouseholderQr().solve(rhs);
      ts.computeStage(&fieldMap);
    }
    ts.computeSolution(&fieldMap);
    double sol = std::accumulate(u.getValues()->begin(), u.getValues()->end(), 0.0)/refEl.getNumNodes();
    if((iter == 0) or (iter == nIters - 1)){
      err2 += std::pow(std::exp(t) - sol, 2)*tstep/2.0;
    } else {
      err2 += std::pow(std::exp(t) - sol, 2)*tstep;
    }
    std::get<2>(*res) = std::sqrt(err2);
  }
};//runRungeKutta

};//rungekutta

TEST_CASE("Testing regression cases for the RungeKutta time schemes", "[regression][RungeKutta]"){
  //list of schemes
  std::vector<RKType> rktypes = {
    //Explicit
    FEuler,
    EMidpoint,
    Heun,
    Kutta3,
    Heun3,
    SSPRK3,
    RK4,
    //Implicit
    BEuler,
    IMidpoint,
    CrankNicolson,
    KS2, //Kraaijevanger and Spijker
    QZ2, //Qin and Zhang
    ALX2, //Alexander 1977
    RK43 //four stage, 3rd order
  };
  std::vector<double> timeSteps = {2e-1, 1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3};
  std::vector<rungekutta::SimRes> results;
  for(int rki = 0; rki < rktypes.size(); rki++){
    for(int dti = 0; dti < timeSteps.size(); dti++){
      results.push_back(std::make_tuple(rktypes[rki], timeSteps[dti], 0.0));
    }
  }
  //std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/RungeKutta/";
  //std::string writeFile = "Breakdown.csv";
  //std::ofstream f;
  //f.open(writePath + writeFile);
  //f << "rktype,timestep,err\n" << std::flush;
  for(int i = 0; i < results.size(); i++){
    rungekutta::runRungeKutta(&results[i]);
    CHECK(std::get<2>(results[i]) < 1.0);
    //RKType rkType = std::get<0>(results[i]);
    //std::string rkstring;
    //if(rkType == BEuler){
      //rkstring = "BEuler";
    //}else if(rkType == ALX2){
      //rkstring = "ALX2";
    //}else if(rkType == IMidpoint){
      //rkstring = "IMidpoint";
    //}else if(rkType == RK4){
      //rkstring = "RK4";
    //}else if(rkType == FEuler){
      //rkstring = "FEuler";
    //}else if(rkType == SSPRK3){
      //rkstring = "SSPRK3";
    //}else if(rkType == EMidpoint){
      //rkstring = "EMidpoint";
    //}else if(rkType == Heun){
      //rkstring = "Heun";
    //}else if(rkType == Kutta3){
      //rkstring = "Kutta3";
    //}else if(rkType == Heun3){
      //rkstring = "Heun3";
    //}else if(rkType == CrankNicolson){
      //rkstring = "CrankNicolson";
    //}else if(rkType == KS2){
      //rkstring = "KS2";
    //}else if(rkType == QZ2){
      //rkstring = "QZ2";
    //}else if(rkType == RK43){
      //rkstring = "RK43";
    //} else{
      //rkstring = "Misc";
    //}
    //f << rkstring << ",";
    //f << std::get<1>(results[i]) << ",";
    //f << std::get<2>(results[i]) << "\n";
    //f << std::flush;
  }
};

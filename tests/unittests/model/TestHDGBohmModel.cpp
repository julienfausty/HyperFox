#include <catch2/catch.hpp>
#include <cmath>
#include "HDGBohmModel.h"

using namespace hfox;

using namespace nGamma;

double tanhTransfer(double x, double y){
  return 0.5*(1.0 - std::tanh((x-y)));
}

std::vector<double> tanhDerivTransfer(double x, double y){
  std::vector<double> res(2, 0);
  double dbuff = 1.0 - std::pow(std::tanh(x-y),2);
  res[0] = -0.5*dbuff;
  res[1] = 0.5*dbuff;
  return res;
}

TEST_CASE("Testing the HDGBohmModel", "[unit][model][HDGBohmModel]"){

  double tol = 1e-5;
  int maxOrd = 5;

  for(int j = 0; j < maxOrd; j++){
    ReferenceElement refEl(1, j+1, "simplex"); 
    SECTION("Testing setup (order = " + std::to_string(j+1) + ")"){
      CHECK_NOTHROW(HDGBohmModel(&refEl));
      HDGBohmModel mod(&refEl);
      CHECK_THROWS(mod.compute());
      CHECK_THROWS(mod.allocate(1));
      CHECK_NOTHROW(mod.allocate(2));
      CHECK_THROWS(mod.compute());
      CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
      std::map<std::string, std::vector<double> > fm;
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Solution"] = std::vector<double>(refEl.getNumNodes()*2, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Trace"] = std::vector<double>(refEl.getNumNodes()*2, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Flux"] = std::vector<double>(refEl.getNumNodes()*4.0, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["ExteriorNormals"] = std::vector<double>(refEl.getNumNodes()*2.0, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["b"] = std::vector<double>(refEl.getNumNodes()*2.0, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["SoundVelocity"] = std::vector<double>(refEl.getNumNodes(), 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Tau"] =std::vector<double>(refEl.getNumNodes()*4.0, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["D"] =std::vector<double>(refEl.getNumNodes()*4.0, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["G"] =std::vector<double>(refEl.getNumNodes()*4.0, 1.0);
      CHECK_NOTHROW(mod.setFieldMap(&fm));
      CHECK_THROWS(mod.compute());
      CHECK_NOTHROW(mod.setTransferFunction(tanhTransfer, tanhDerivTransfer));
    };

    SECTION("Test Dirichlet condition (order=" + std::to_string(j+1) + ")"){
      HDGBohmModel mod(&refEl);
      int nNodes = refEl.getNumNodes();
      std::vector< std::vector<double> > elementNodes(nNodes, std::vector<double>(2, 0));
      for(int k = 0; k < nNodes; k++){
        elementNodes[k][0] = refEl.getNodes()->at(k)[0];
      }
      CHECK_NOTHROW(mod.setElementNodes(&elementNodes));
      mod.allocate(2);
      std::map<std::string, std::vector<double> > fm;
      fm["Solution"] = std::vector<double>(nNodes*2, 1.0);
      fm["Trace"] = std::vector<double>(nNodes*2, 1.0);
      fm["Flux"] = std::vector<double>(nNodes*4, 0.0);
      fm["Tau"] = std::vector<double>(nNodes*4, 0.0);
      fm["D"] = std::vector<double>(nNodes*4, 0.0);
      fm["G"] = std::vector<double>(nNodes*4, 0.0);
      fm["ExteriorNormals"] = std::vector<double>(nNodes*2.0, 0.0);
      for(int k = 0; k < nNodes; k++){
        fm["ExteriorNormals"][2*k+1] = -1.0;
        EMap<EMatrix>(fm["Tau"].data() + 4*k, 2, 2) = EMatrix::Identity(2,2);
        EMap<EMatrix>(fm["D"].data() + 4*k, 2, 2) = EMatrix::Identity(2,2);
        EMap<EMatrix>(fm["G"].data() + 4*k, 2, 2) = EMatrix::Identity(2,2);
      }
      fm["b"] = fm["ExteriorNormals"];
      fm["SoundVelocity"] = std::vector<double>(nNodes, 10.0);
      mod.setFieldMap(&fm);
      mod.setTransferFunction(tanhTransfer, tanhDerivTransfer);
      CHECK_NOTHROW(mod.compute());
      EVector testSolution(8*nNodes);
      for(int k = 0; k < nNodes; k++){
        testSolution(2*k) = 1.0;
        testSolution(2*k+1) = 10.0;
        testSolution(2*nNodes + 4*k) = 1.0;
        testSolution(2*nNodes + 4*k + 1) = 0.0;
        testSolution(2*nNodes + 4*k + 2) = 0.0;
        testSolution(2*nNodes + 4*k + 3) = 1.0;
        testSolution(6*nNodes + 2*k) = 1.0;
        testSolution(6*nNodes + 2*k+1) = 10.0;
      }
      CHECK((*(mod.getLocalMatrix())*testSolution - *(mod.getLocalRHS())).norm() < tol);
    };

    SECTION("Test Neumann condition (order=" + std::to_string(j+1) + ")"){
      HDGBohmModel mod(&refEl);
      int nNodes = refEl.getNumNodes();
      std::vector< std::vector<double> > elementNodes(nNodes, std::vector<double>(2, 0));
      for(int k = 0; k < nNodes; k++){
        elementNodes[k][0] = refEl.getNodes()->at(k)[0];
      }
      CHECK_NOTHROW(mod.setElementNodes(&elementNodes));
      mod.allocate(2);
      std::map<std::string, std::vector<double> > fm;
      fm["Solution"] = std::vector<double>(nNodes*2, 1.0);
      fm["Trace"] = std::vector<double>(nNodes*2, 1.0);
      fm["Flux"] = std::vector<double>(nNodes*4, 0.0);
      fm["ExteriorNormals"] = std::vector<double>(nNodes*2.0, 0.0);
      fm["Tau"] = std::vector<double>(nNodes*4, 0.0);
      fm["D"] = std::vector<double>(nNodes*4, 0.0);
      fm["G"] = std::vector<double>(nNodes*4, 0.0);
      for(int k = 0; k < nNodes; k++){
        fm["Solution"][2*k+1] = 10.0;
        fm["Trace"][2*k+1] = 10.0;
        fm["ExteriorNormals"][2*k+1] = -1.0;
        EMap<EMatrix>(fm["Tau"].data() + 4*k, 2, 2) = EMatrix::Identity(2,2);
        EMap<EMatrix>(fm["D"].data() + 4*k, 2, 2) = EMatrix::Identity(2,2);
        EMap<EMatrix>(fm["G"].data() + 4*k, 2, 2) = EMatrix::Identity(2,2);
      }
      fm["b"] = fm["ExteriorNormals"];
      fm["SoundVelocity"] = std::vector<double>(nNodes, 1.0);
      mod.setFieldMap(&fm);
      mod.setTransferFunction(tanhTransfer, tanhDerivTransfer);
      CHECK_NOTHROW(mod.compute());
      EVector testSolution(8*nNodes);
      for(int k = 0; k < nNodes; k++){
        testSolution(2*k) = 1.0;
        testSolution(2*k+1) = 1.0;
        testSolution(2*nNodes + 4*k) = 1.0;
        testSolution(2*nNodes + 4*k + 1) = 1.0;
        testSolution(2*nNodes + 4*k + 2) = 0.0;
        testSolution(2*nNodes + 4*k + 3) = 0.0;
        testSolution(6*nNodes + 2*k) = 1.0;
        testSolution(6*nNodes + 2*k+1) = 1.0;
      }
      CHECK((*(mod.getLocalMatrix())*testSolution - *(mod.getLocalRHS())).norm() < tol);
    };

  }

};

#include <catch2/catch.hpp>
#include <string>
#include <map>
#include <cmath>
#include "Operator.h"
#include "HDGBase.h"
#include "HDGDiffusion.h"
#include "HDGLaplaceModel.h"

using namespace hfox;

TEST_CASE("Testing the HDGLaplaceModel", "[unit][model][HDGLaplaceModel]"){

  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      SECTION("Testing setup (dim="+std::to_string(i+1)+", order="+std::to_string(j+1)+")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        HDGLaplaceModel hdgLapMod(&refEl);
        CHECK_THROWS(hdgLapMod.compute());
        std::map<std::string, std::vector<double> > fm;
        CHECK_THROWS(hdgLapMod.setFieldMap(&fm));
        std::vector<double> taus(refEl.getNumFaces()*(refEl.getFaceElement()->getNumNodes()), 1.0);
        fm["Tau"] = taus;
        CHECK_NOTHROW(hdgLapMod.setFieldMap(&fm));
        CHECK_THROWS(hdgLapMod.compute());
        CHECK_NOTHROW(hdgLapMod.setElementNodes(refEl.getNodes()));
        CHECK_THROWS(hdgLapMod.compute());
        CHECK_NOTHROW(hdgLapMod.allocate(1));
        CHECK_NOTHROW(hdgLapMod.compute());
      };

      ReferenceElement refEl(i+1, j+1, "simplex");
      HDGLaplaceModel hdgLapMod(&refEl);
      std::map<std::string, std::vector<double> > fm;
      std::vector<double> taus(refEl.getNumFaces()*(refEl.getFaceElement()->getNumNodes()), 1.0);
      fm["Tau"] = taus;
      hdgLapMod.setFieldMap(&fm);
      hdgLapMod.setElementNodes(refEl.getNodes());
      hdgLapMod.allocate(1);
      hdgLapMod.compute();
      double tol = 1e-12;

      SECTION("Testing Matrix (dim="+std::to_string(i+1)+", order="+std::to_string(j+1)+")"){
        HDGBase baseOp(&refEl);
        HDGDiffusion diffOp(&refEl);
        baseOp.allocate(1);
        diffOp.allocate(1);
        baseOp.setTau(taus);
        std::vector<EMatrix> elJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
        std::vector<EMatrix> jacs(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs());
        std::vector<EMatrix> invJacs(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs());
        std::vector<double> dV(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs(), 0.0);
        std::copy(elJacs.begin(), elJacs.end(), jacs.begin());
        std::copy(elJacs.begin(), elJacs.end(), invJacs.begin());
        std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
        std::vector<EMatrix> faceJacs(refEl.getFaceElement()->getNumIPs(), EMatrix::Zero(i, i+1));
        std::vector<double> facedV(refEl.getFaceElement()->getNumIPs(), 0.0);
        std::vector< std::vector<double> > faceNodes(refEl.getFaceElement()->getNumNodes(), std::vector<double>(i+1, 0.0));
        for(int k = 0; k < refEl.getNumFaces(); k++){
          for(int l = 0; l < refEl.getFaceElement()->getNumNodes(); l++){
            faceNodes[l] = refEl.getNodes()->at(refEl.getFaceNodes()->at(k)[l]);
          }
          int offset = refEl.getNumIPs() + k*(refEl.getFaceElement()->getNumIPs());
          faceJacs = Operator::calcJacobians(faceNodes, refEl.getFaceElement());
          std::copy(faceJacs.begin(), faceJacs.end(), jacs.begin() + offset);
          facedV = Operator::calcMeasure(Operator::calcDetJacobians(faceJacs), refEl.getFaceElement());
          std::copy(facedV.begin(), facedV.end(), dV.begin() + offset);
          faceJacs = Operator::calcInvJacobians(faceJacs);
          std::copy(faceJacs.begin(), faceJacs.end(), invJacs.begin() + offset);
        }
        baseOp.calcNormals(*(refEl.getNodes()), jacs);
        baseOp.assemble(dV, invJacs);
        diffOp.setFromBase(baseOp.getNormals());
        diffOp.assemble(dV, invJacs);
        EMatrix testMat = *(baseOp.getMatrix()) + *(diffOp.getMatrix());
        testMat -= *(hdgLapMod.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
      };
      SECTION("Testing RHS (dim="+std::to_string(i+1)+", order="+std::to_string(j+1)+")"){
        CHECK(std::abs((hdgLapMod.getLocalRHS())->sum()) < tol);
      }
    }
  }
};

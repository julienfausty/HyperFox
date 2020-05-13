#include <catch2/catch.hpp>
#include "HDGDiffusionSource.h"
#include "Euler.h"

using namespace hfox;

double testSrc(const std::vector<double> & v){
  return v[0];
};

TEST_CASE("Testing the HDGDiffusionSource model", "[unit][model][HDGDiffusionSource]"){

  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      SECTION("Testing setup (dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(HDGDiffusionSource(&refEl));
        HDGDiffusionSource mod(&refEl);
        CHECK_THROWS(mod.compute());
        std::map<std::string, std::vector<double> > fm;
        CHECK_THROWS(mod.setFieldMap(&fm));
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        std::vector<double> buff(refEl.getNumNodes(), 0.0);
        fm["Tau"] = buff;
        CHECK_NOTHROW(mod.setFieldMap(&fm));
        CHECK_THROWS(mod.compute());
        CHECK_THROWS(mod.setSourceFunction(testSrc));
        CHECK_NOTHROW(mod.allocate(1));
        CHECK_NOTHROW(mod.setSourceFunction(testSrc));
      };

      HDGDiffusionSource mod(&refEl);
      std::map<std::string, std::vector<double> > fm;
      std::vector<double> tau(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()), 1.0);
      std::vector<double> sol(refEl.getNumNodes(), 0.0);
      std::iota(sol.begin(), sol.end(), 0.0);
      fm["Tau"] = tau;
      fm["Solution"] = sol;
      mod.setFieldMap(&fm);
      mod.allocate(1);
      mod.setSourceFunction(testSrc);
      HDGBase baseOp(&refEl);
      Source srcOp(&refEl);
      baseOp.allocate(1);
      baseOp.setTau(tau);
      srcOp.allocate(1);
      srcOp.setSourceFunction(testSrc);
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
      srcOp.calcSource(*(refEl.getNodes()));
      srcOp.assemble(dV, invJacs);
      int startU = 0;
      int startQ = refEl.getNumNodes();
      int startL = refEl.getNumNodes()*(i+2);
      int lenU = refEl.getNumNodes();
      int lenQ = refEl.getNumNodes() * (i+1);
      int lenL = (refEl.getFaceElement()->getNumNodes())*(refEl.getNumFaces());
      //Suq
      Convection convOp(&refEl);
      convOp.allocate(1);
      EMatrix Suq(lenU, lenQ);
      for(int k = 0; k < i+1; k++){
        EVector unitBase = EVector::Zero(i+1);
        unitBase[k] = 1.0;
        std::vector<EVector> velocity(lenU, unitBase);
        convOp.setVelocity(velocity);
        convOp.assemble(*(refEl.getIPWeights()), invJacs);
        for(int l = 0 ; l < lenU; l++){
          Suq.col(l*(i+1) + k) = (convOp.getMatrix()->row(l)).transpose();
        }
      }
      for(int iFace = 0; iFace < refEl.getNumFaces(); iFace++){
        for(int k = 0; k < refEl.getFaceElement()->getNumNodes(); k++){
          int rowInd = (refEl.getFaceNodes()->at(iFace))[k];
          Suq.row(rowInd) -= baseOp.getMatrix()->block(startL + iFace*(refEl.getFaceElement()->getNumNodes()) + k, startQ, 1, lenQ);
        }
      }
      SECTION("Poisson equation in reference element"){
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod.compute());
        EMatrix testMat = *(baseOp.getMatrix());
        testMat.block(startU, startQ, lenU, lenQ) += Suq;
        testMat -= *(mod.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = (*(mod.getLocalRHS()));
        testVec.segment(0, lenU) -= ((EVector) *(srcOp.getMatrix()));
        CHECK(testVec.dot(testVec) < tol);
      };
      HDGDiffusionSource mod2(&refEl);
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
      mod2.setSourceFunction(testSrc);
      SECTION("Diffusion source equation in reference element"){
        CHECK_NOTHROW(mod2.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod2.compute());
        EMatrix testMat = *(baseOp.getMatrix());
        testMat.block(startU, startQ, lenU, lenQ) += 3.0*Suq;
        testMat.block(0, 0, lenU, testMat.cols()) *= timeStep;
        testMat.block(0, 0, lenU, lenU) += (*(mass.getMatrix()));
        testMat -= *(mod2.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = timeStep*((EVector) *(srcOp.getMatrix()));
        testVec += (*(mass.getMatrix()))*EMap<EVector>(sol.data(), sol.size());
        testVec -= (*(mod2.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };
    }
  }
};

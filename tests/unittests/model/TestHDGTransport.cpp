#include <catch2/catch.hpp>
#include "HDGTransport.h"
#include "Euler.h"

using namespace hfox;

TEST_CASE("Testing the HDGTransport model", "[unit][model][HDGTransport]"){
  
  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      SECTION("Testing setup (dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(HDGTransport(&refEl));
        HDGTransport mod(&refEl);
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.allocate(1));
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        std::map<std::string, std::vector<double> > fm;
        CHECK_THROWS(mod.setFieldMap(&fm));
        fm["Velocity"] = std::vector<double>(refEl.getNumNodes()*(i+1), 1.0);
        CHECK_THROWS(mod.setFieldMap(&fm));
        fm["Tau"] = std::vector<double>(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()), 1.0);
        CHECK_NOTHROW(mod.setFieldMap(&fm));
      };

      HDGTransport mod(&refEl);
      std::map<std::string, std::vector<double> > fm;
      std::vector<double> tau(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()), 1.0);
      std::vector<double> sol(refEl.getNumNodes(), 0.0);
      std::vector<double> vel(refEl.getNumNodes()*(i+1), 2.0);
      std::iota(sol.begin(), sol.end(), 0.0);
      fm["Tau"] = tau;
      fm["Solution"] = sol;
      fm["Velocity"] = vel;
      mod.setFieldMap(&fm);
      mod.allocate(1);
      HDGBase baseOp(&refEl);
      std::vector<EVector> velVectors(refEl.getNumNodes(), EVector::Constant(i+1, 2.0));
      HDGConvection convOp(&refEl);
      baseOp.allocate(1);
      baseOp.setTau(tau);
      convOp.allocate(1);
      convOp.setVelocity(velVectors);
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
      convOp.setFromBase(baseOp.getNormals());
      convOp.assemble(dV, invJacs);
      SECTION("Grad ortho equation in reference element(dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod.compute());
        EMatrix testMat = *(convOp.getMatrix()) + (*(baseOp.getMatrix()));
        testMat -= *(mod.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = (*(mod.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };

      HDGTransport mod2(&refEl);
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
        EMatrix testMat = *(convOp.getMatrix()) + (*(baseOp.getMatrix()));
        testMat.block(0, 0, refEl.getNumNodes(), testMat.cols()) *= timeStep;
        testMat.block(0, 0, refEl.getNumNodes(), refEl.getNumNodes()) += (*(mass.getMatrix()));
        testMat -= *(mod2.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = EVector::Zero(testMat.rows());
        testVec.segment(0, refEl.getNumNodes()) = (*(mass.getMatrix()))*EMap<EVector>(sol.data(), sol.size());
        testVec -= (*(mod2.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };
    }
  }


};

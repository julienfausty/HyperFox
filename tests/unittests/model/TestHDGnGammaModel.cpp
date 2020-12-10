#include <catch2/catch.hpp>
#include "HDGnGammaModel.h"
#include "HDGConvection.h"
#include "HDGDiffusion.h"
#include "Reaction.h"
#include "Source.h"
#include "Euler.h"

using namespace hfox;
using namespace nGamma;

TEST_CASE("Testing the HDGnGammaModel", "[unit][model][HDGnGammaModel]"){

  double tol = 1e-12;
  int maxOrd = 7;

  for(int j = 0; j < maxOrd; j++){

    ReferenceElement refEl(2, j+1, "simplex");
    SECTION("Testing setup (order=" + std::to_string(j+1) + ")"){
      CHECK_NOTHROW(HDGnGammaModel(&refEl));
      HDGnGammaModel mod(&refEl);
      CHECK_THROWS(mod.compute());
      CHECK_THROWS(mod.allocate(3));
      CHECK_NOTHROW(mod.allocate(2));
      CHECK_THROWS(mod.compute());
      CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
      std::map<std::string, std::vector<double> > fm;
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["BufferSolution"] = std::vector<double>(refEl.getNumNodes()*2, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Tau"] = std::vector<double>(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * std::pow(2, 2), 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      //fm["Trace"] = std::vector<double>(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * 2, 1.0);
      //CHECK_THROWS(mod.setFieldMap(&fm));
      fm["D"] = std::vector<double>(refEl.getNumNodes()*4, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["G"] = std::vector<double>(refEl.getNumNodes()*4, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["b"] = std::vector<double>(refEl.getNumNodes()*3, 1.0);
      CHECK_NOTHROW(mod.setFieldMap(&fm));
    };

    //setup jacobians etc.
    std::vector<EMatrix> elJacs(refEl.getNumIPs(), EMatrix::Identity(2, 2));
    std::vector<EMatrix> jacs(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs());
    std::vector<EMatrix> invJacs(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs());
    std::vector<double> dV(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs(), 0.0);
    std::copy(elJacs.begin(), elJacs.end(), jacs.begin());
    std::copy(elJacs.begin(), elJacs.end(), invJacs.begin());
    std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
    std::vector<EMatrix> faceJacs(refEl.getFaceElement()->getNumIPs(), EMatrix::Zero(1, 2));
    std::vector<double> facedV(refEl.getFaceElement()->getNumIPs(), 0.0);
    std::vector< std::vector<double> > faceNodes(refEl.getFaceElement()->getNumNodes(), std::vector<double>(2, 0.0));
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

    //setup field map
    std::map<std::string, std::vector<double> > fm;
    std::vector<double> tau(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * std::pow(2, 2), 0.0);
    std::vector<double> trace(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * 2, 0.0);
    std::vector<double> sol(refEl.getNumNodes() * 2, 2.0);
    std::vector<double> diff(refEl.getNumNodes()*4, 0);
    std::vector<double> b(refEl.getNumNodes() * 2, 1.0);
    for(int k = 0; k < refEl.getNumNodes(); k++){
      EMap<EMatrix>(diff.data() + 4*k, 2, 2) = EMatrix::Identity(2, 2)*3.0;
    }
    for(int k = 0; k < refEl.getNumFaces()*refEl.getFaceElement()->getNumNodes(); k++){
      EMap<EMatrix>(tau.data() + 4*k, 2, 2) = EMatrix::Identity(2, 2);
    }
    int nFaces = refEl.getNumFaces();
    int nNodesFc = refEl.getFaceElement()->getNumNodes();
    for(int iFace = 0; iFace < nFaces; iFace++){
      for(int iN = 0; iN < nNodesFc; iN++){
        for(int d = 0; d < 2; d++){
          trace[(iFace*nNodesFc + iN)*2 + d] = 2.0;
        }
      }
    }
    fm["Tau"] = tau;
    fm["BufferSolution"] = sol;
    //fm["Trace"] = trace;
    fm["D"] = diff;
    fm["G"] = diff;
    fm["b"] = b;

    //parameters
    nGammaParams paramSet;
    paramSet.soundSpeed = 0.5;

    //setup operators
    HDGBase baseOp(&refEl);
    std::vector<HDGConvection*> convectionOps(4);
    std::vector<HDGDiffusion*> diffusionOps(2);
    std::vector< std::vector<EVector> > allVelocities(4, std::vector<EVector>(refEl.getNumNodes(), EVector::Zero(2)));
    std::vector< std::vector<EMatrix> > allDiffs(2, std::vector<EMatrix>(refEl.getNumNodes(), EMatrix::Identity(2, 2)*3.0));
    for(int k = 0; k < refEl.getNumNodes(); k++){
      allVelocities[1][k] = EVector::Constant(2, 1.0);//b
      allVelocities[2][k] = (std::pow(paramSet.soundSpeed, 2) - 1.0)*EVector::Constant(2, 1.0);//(cs^2 - U1^2/U0^2)b
      allVelocities[3][k] = 2.0*EVector::Constant(2, 1.0);//(2U1/U0)b
    }
    baseOp.allocate(2);
    baseOp.setTau(tau);
    baseOp.calcNormals(*(refEl.getNodes()), jacs);
    baseOp.assemble(dV, invJacs);

    for(int k = 0; k < 4; k++){
      convectionOps[k] = new HDGConvection(&refEl);
      convectionOps[k]->allocate(1);
      convectionOps[k]->setFromBase(baseOp.getNormals());
      convectionOps[k]->setVelocity(allVelocities[k]);
      convectionOps[k]->assemble(dV, invJacs);
    }
    for(int k = 0; k < 2; k++){
      diffusionOps[k] = new HDGDiffusion(&refEl);
      diffusionOps[k]->allocate(1);
      diffusionOps[k]->setFromBase(baseOp.getNormals());
      diffusionOps[k]->setDiffusionTensor(allDiffs[k]);
      diffusionOps[k]->assemble(dV, invJacs);
    }

    //assemble stiffness
    int stiffSize = sol.size() * 3 + trace.size();
    EMatrix stiffness(stiffSize, stiffSize);
    stiffness = *(baseOp.getMatrix());
    for(int k = 0; k < stiffSize/2; k++){
      for(int l = 0; l < stiffSize/2; l++){
        stiffness(k*2, l*2) += (*(convectionOps[0]->getMatrix()))(k, l) + (*(diffusionOps[0]->getMatrix()))(k, l);
        stiffness(k*2, l*2+1) += (*(convectionOps[1]->getMatrix()))(k, l);
        stiffness(k*2 + 1, l*2) += (*(convectionOps[2]->getMatrix()))(k, l);
        stiffness(k*2 + 1, l*2 + 1) += (*(convectionOps[3]->getMatrix()))(k, l) + (*(diffusionOps[1]->getMatrix()))(k, l);
      }
    }

    SECTION("Testing stationary assembly (order = " + std::to_string(j+1) + ")"){
      HDGnGammaModel mod(&refEl);
      mod.setParams(paramSet);
      mod.allocate(2);
      mod.setFieldMap(&fm);
      mod.setElementNodes(refEl.getNodes());
      CHECK_NOTHROW(mod.compute());
      EMatrix testMat = stiffness;
      testMat -= *(mod.getLocalMatrix());
      //std::cout << "test calculated:" << std::endl;
      //std::cout << stiffness << std::endl;
      //std::cout << "model calculated:" << std::endl;
      //std::cout << *(mod.getLocalMatrix()) << std::endl;
      //std::cout << "subtraction:" << std::endl;
      //std::cout << testMat << std::endl;
      CHECK((testMat.transpose() * testMat).sum() < tol);
      EVector testVec = *(mod.getLocalRHS());
      CHECK(testVec.dot(testVec) < tol);
    };

    for(int k = 0; k < convectionOps.size(); k++){
      delete convectionOps[k];
    }

    for(int k = 0; k < diffusionOps.size(); k++){
      delete diffusionOps[k];
    }
  }

};

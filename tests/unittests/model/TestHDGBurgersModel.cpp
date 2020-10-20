#include <catch2/catch.hpp>
#include "HDGBurgersModel.h"
#include "Euler.h"

using namespace hfox;

TEST_CASE("Testing the HDGBurgersModel", "[unit][model][HDGBurgersModel]"){
  
  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      SECTION("Testing setup (dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(HDGBurgersModel(&refEl));
        HDGBurgersModel mod(&refEl);
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.allocate(i+1));
        CHECK_THROWS(mod.compute());
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        std::map<std::string, std::vector<double> > fm;
        CHECK_THROWS(mod.setFieldMap(&fm));
        fm["Solution"] = std::vector<double>(refEl.getNumNodes()*(i+1), 1.0);
        CHECK_THROWS(mod.setFieldMap(&fm));
        fm["Tau"] = std::vector<double>(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * std::pow(i+1, 2), 1.0);
        CHECK_THROWS(mod.setFieldMap(&fm));
        fm["Trace"] = std::vector<double>(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * (i+1), 1.0);
        CHECK_NOTHROW(mod.setFieldMap(&fm));
      };

      HDGBurgersModel mod(&refEl);
      std::map<std::string, std::vector<double> > fm;
      std::vector<double> tau(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * std::pow(i+1, 2), 0.0);
      std::vector<double> trace(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * (i+1), 0.0);
      std::vector<double> sol(refEl.getNumNodes() * (i+1), 0.0);
      std::vector<double> buffSol(refEl.getNumNodes() * (i+1), 2.0);
      std::vector<double> diff(refEl.getNumNodes()*(i+1)*(i+1), 0);
      for(int k = 0; k < refEl.getNumNodes(); k++){
        EMap<EMatrix>(diff.data() + (i+1)*(i+1)*k, (i+1), (i+1)) = EMatrix::Identity((i+1), (i+1))*3.0;
      }
      for(int k = 0; k < refEl.getNumFaces()*refEl.getFaceElement()->getNumNodes(); k++){
        EMap<EMatrix>(tau.data() + (i+1)*(i+1)*k, (i+1), (i+1)) = EMatrix::Identity((i+1), (i+1));
      }
      std::fill(sol.begin(), sol.end(), 2.0);
      int nFaces = refEl.getNumFaces();
      int nNodesFc = refEl.getFaceElement()->getNumNodes();
      for(int iFace = 0; iFace < nFaces; iFace++){
        for(int iN = 0; iN < nNodesFc; iN++){
          for(int d = 0; d < (i+1); d++){
            trace[(iFace*nNodesFc + iN)*(i+1) + d] = 2.0;
          }
        }
      }
      fm["Tau"] = tau;
      fm["Solution"] = sol;
      fm["DiffusionTensor"] = diff;
      fm["Trace"] = trace;
      mod.setFieldMap(&fm);
      mod.allocate(i+1);
      HDGBase baseOp(&refEl);
      std::vector<EVector> velVectors(refEl.getNumNodes(), EVector::Constant(i+1, 2.0));
      std::vector<EVector> traceVectors(nFaces*nNodesFc, EVector::Constant(i+1, 2.0));
      HDGUNabU convOp(&refEl);
      std::vector<EMatrix> diffTensors(refEl.getNumNodes(), EMatrix::Identity(i+1, i+1)*3.0);
      HDGDiffusion diffOp(&refEl);
      baseOp.allocate(i+1);
      baseOp.setTau(tau);
      convOp.allocate(i+1);
      convOp.setSolution(velVectors);
      convOp.setTrace(traceVectors);
      diffOp.allocate(i+1);
      diffOp.setDiffusionTensor(diffTensors);
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
      diffOp.setFromBase(baseOp.getNormals());
      diffOp.assemble(dV, invJacs);
      SECTION("No time operator in reference element(dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(mod.compute());
        EMatrix testMat = (*(baseOp.getMatrix()));
        testMat += *(convOp.getMatrix()) + *(diffOp.getMatrix());
        testMat -= *(mod.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = (*(mod.getLocalRHS()));
        testVec -= (*(convOp.getRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };

      HDGBurgersModel mod2(&refEl);
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
        EMatrix testMat = *(convOp.getMatrix()) + (*(baseOp.getMatrix())) + *(diffOp.getMatrix());
        testMat.block(0, 0, refEl.getNumNodes()*(i+1), testMat.cols()) *= timeStep;
        testMat.block(0, 0, refEl.getNumNodes()*(i+1), refEl.getNumNodes()*(i+1)) += (*(mass.getMatrix()));
        testMat -= *(mod2.getLocalMatrix());
        CHECK((testMat.transpose() * testMat).sum() < tol);
        EVector testVec = (*(convOp.getRHS()));
        testVec.segment(0, refEl.getNumNodes()*(i+1)) = (*(mass.getMatrix()))*EMap<EVector>(sol.data(), sol.size()); 
        testVec -= (*(mod2.getLocalRHS()));
        CHECK(testVec.dot(testVec) < tol);
      };
    }
  }


};

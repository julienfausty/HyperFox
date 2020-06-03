#include <catch2/catch.hpp>
#include "HDGConvection.h"
#include "Mass.h"
#include "TestUtils.h"

using namespace hfox;

TEST_CASE("Testing the HDGConvection operator", "[unit][operator][HDGConvection]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(2, 1, "simplex");
    HDGConvection op(&refEl);
    int n = refEl.getNumIPs()+refEl.getNumFaces()*(refEl.getFaceElement()->getNumIPs());
    std::vector<double> dummydV(n, 0.0);
    std::vector<EMatrix> dummyJac(n, EMatrix::Identity(2,2));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(op.allocate(1));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
  };
  
  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      HDGConvection op(&refEl);
      op.allocate(1);
      std::vector< std::vector<double> > refNorms = TestUtils::getRefNormals(i+1);
      int nNodesEl = refEl.getNumNodes();
      int nIPsEl = refEl.getNumIPs();
      int nFaces = refEl.getNumFaces();
      const ReferenceElement * fEl = refEl.getFaceElement();
      int nNodesPFc = fEl->getNumNodes();
      int nIPsPFc = fEl->getNumIPs();
      SECTION("Testing assembly in reference element (dim=" + std::to_string(i+1) + ", order=" + std::to_string(j+1) + ")"){
        std::vector<double> dV(nIPsEl + nFaces*nIPsPFc, 0.0);
        std::vector<EMatrix> jacs(nIPsEl + nFaces*nIPsPFc);
        std::fill(jacs.begin(), jacs.begin() + nIPsEl, EMatrix::Identity(i+1, i+1));
        std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
        std::vector<EVector> normals(nFaces*nIPsPFc, EVector::Zero(i+1));
        std::vector< std::vector<double> > faceNodes(nNodesPFc, std::vector<double>(i+1));
        std::vector<EMatrix> fcMats(nIPsPFc);
        std::vector<double> fcdV(nIPsPFc, 0.0);
        Mass fcMass(refEl.getFaceElement());
        fcMass.allocate(1);
        std::vector<EMatrix> fcMasses(nFaces, EMatrix::Zero(nNodesPFc, nNodesPFc));
        for(int k = 0; k < nFaces; k++){
          std::fill(normals.begin() + k*nIPsPFc, normals.begin() + (k+1)*nIPsPFc, EMap<EVector>(refNorms[k].data(), refNorms[k].size()));
          for(int l = 0; l < nNodesPFc; l++){
            faceNodes[l] = refEl.getNodes()->at(refEl.getFaceNodes()->at(k)[l]);
          }
          fcMats = Operator::calcJacobians(faceNodes, fEl);
          fcdV = Operator::calcMeasure(Operator::calcDetJacobians(fcMats), fEl);
          fcMats = Operator::calcInvJacobians(fcMats);
          std::copy(fcdV.begin(), fcdV.end(), dV.begin()+ nIPsEl + k*nIPsPFc);
          std::copy(fcMats.begin(), fcMats.end(), jacs.begin() + nIPsEl + k*nIPsPFc);
          fcMass.assemble(fcdV, fcMats);
          fcMasses[k] = *(fcMass.getMatrix());
        }
        std::vector<EVector> vels(nNodesEl, EVector::Constant(i+1, 2.0));
        CHECK_NOTHROW(op.setFromBase(&normals));
        CHECK_NOTHROW(op.setVelocity(vels));
        op.setFromBase(&normals);
        CHECK_NOTHROW(op.assemble(dV, jacs));
        Convection convOp(&refEl);
        convOp.allocate(1);
        convOp.setVelocity(vels);
        convOp.assemble(dV, jacs);
        EMatrix Suu = EMatrix::Zero(nNodesEl, nNodesEl);
        EMatrix Sll = EMatrix::Zero(nNodesPFc*nFaces, nNodesPFc*nFaces);
        EMatrix Sul = EMatrix::Zero(nNodesEl, nNodesPFc*nFaces);
        Suu = -((convOp.getMatrix())->transpose());
        for(int iFace = 0; iFace < nFaces; iFace++){
          double prefactor = EVector::Constant(i+1, 2.0).dot(EMap<const EVector>(refNorms[iFace].data(), refNorms[iFace].size()));
          for(int k = 0; k < nNodesPFc; k++){
            int rowInd = (refEl.getFaceNodes()->at(iFace))[k];
            for(int l = 0; l < nNodesPFc; l++){
              int colInd = (refEl.getFaceNodes()->at(iFace))[l];
              Sul(rowInd, iFace*nNodesPFc + l) += fcMasses[iFace](k, l)*prefactor;
              Sll(iFace*nNodesPFc + k, iFace*nNodesPFc + l) += fcMasses[iFace](k, l)*prefactor;
            }
          }
        }
        EMatrix testMat = *(op.getMatrix());
        testMat.block(0, 0, nNodesEl, nNodesEl) -= Suu;
        testMat.block(0, nNodesEl*(i+2), nNodesEl, nNodesPFc*nFaces) -= Sul;
        testMat.block(nNodesEl*(i+2), nNodesEl*(i+2), nNodesPFc*nFaces, nNodesPFc*nFaces) -= Sll;
        CHECK((testMat.transpose() * testMat).sum() < tol);
      };
    }
  }

};

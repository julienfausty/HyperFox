#include <catch2/catch.hpp>
#include "HDGDiffusion.h"
#include "Convection.h"
#include "Mass.h"
#include "TestUtils.h"

using namespace hfox;

TEST_CASE("Testing HDGDiffusion operator", "[unit][operator][HDGDiffusion]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(2, 1, "simplex");
    HDGDiffusion diff(&refEl);
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(diff.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(diff.allocate(1));
    CHECK_THROWS(diff.assemble(dummydV, dummyJac));
  };

  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      HDGDiffusion diff(&refEl);
      diff.allocate(1);
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
        CHECK_NOTHROW(diff.setFromBase(&normals));
        CHECK_NOTHROW(diff.assemble(dV, jacs));
        Convection convOp(&refEl);
        convOp.allocate(1);
        EMatrix Suq = EMatrix::Zero(nNodesEl, nNodesEl*(i+1));
        for(int k = 0; k < (i+1); k++){
          EVector unitBase = EVector::Zero(i+1);
          unitBase[k] = 1.0;
          std::vector<EVector> velocity(nNodesEl, unitBase);
          convOp.setVelocity(velocity);
          convOp.assemble(*(refEl.getIPWeights()), jacs);
          for(int l = 0 ; l < nNodesEl; l++){
            Suq.col(l*(i+1) + k) += (convOp.getMatrix()->row(l)).transpose();
          }
        }
        for(int iFace = 0; iFace < nFaces; iFace++){
          for(int k = 0; k < nNodesPFc; k++){
            int rowInd = (refEl.getFaceNodes()->at(iFace))[k];
            for(int l = 0; l < nNodesPFc; l++){
              int colInd = (refEl.getFaceNodes()->at(iFace))[l];
              for(int d = 0; d < i+1; d++){
                Suq(rowInd, colInd*(i+1) + d) -= fcMasses[iFace](k, l)*refNorms[iFace][d];
              }
            }
          }
        }
        EMatrix testMat = diff.getMatrix()->block(0, nNodesEl, nNodesEl, nNodesEl*(i+1)) - Suq;
        CHECK((testMat.transpose() * testMat).sum() < tol);
        std::vector<EMatrix> diffTensor(nNodesEl, 3.0*EMatrix::Identity(i+1, i+1)); 
        CHECK_NOTHROW(diff.setDiffusionTensor(diffTensor));
        diff.assemble(dV, jacs);
        testMat = diff.getMatrix()->block(0, nNodesEl, nNodesEl, nNodesEl*(i+1)) - 3.0*Suq;
        CHECK((testMat.transpose() * testMat).sum() < tol);
      };
    }
  }

};

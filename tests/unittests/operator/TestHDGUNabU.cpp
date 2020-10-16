#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include <algorithm>
#include "TestUtils.h"
#include "HDGUNabU.h"

using namespace hfox;

TEST_CASE("Testing HDGUNabU operator", "[unit][operator][nonlinear][HDGUNabU]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    HDGUNabU op(&refEl);
    std::vector<EVector> velocity(refEl.getNumNodes(), EVector::Zero(1));
    std::vector<EVector> trace(refEl.getFaceElement()->getNumNodes()*refEl.getNumFaces(), EVector::Zero(1));
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(op.allocate(1));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(op.setSolution(velocity));
    CHECK_NOTHROW(op.setTrace(trace));
  };

  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      for(int n = 0; n < 2*j/3; n++){
        SECTION("Testing HDGUNabU assembly in reference element: (dim = " + std::to_string(i+1) + ", order = " + std::to_string(j+1) + ", n= " + std::to_string(n+1) + ")"){
          if((n+1) % 2 != 0){
            ReferenceElement refEl(i+1, j+1, "simplex");
            std::vector< std::vector<double> > refNorms = TestUtils::getRefNormals(i+1);
            int nNodesEl = refEl.getNumNodes();
            int nIPsEl = refEl.getNumIPs();
            int nFaces = refEl.getNumFaces();
            const ReferenceElement * fEl = refEl.getFaceElement();
            int nNodesPFc = fEl->getNumNodes();
            int nIPsPFc = fEl->getNumIPs();
            std::vector<double> dV(nIPsEl + nFaces*nIPsPFc, 0.0);
            std::vector<EMatrix> jacs(nIPsEl + nFaces*nIPsPFc);
            std::fill(jacs.begin(), jacs.begin() + nIPsEl, EMatrix::Identity(i+1, i+1));
            std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
            std::vector<EVector> normals(nFaces*nIPsPFc, EVector::Zero(i+1));
            std::vector< std::vector<double> > faceNodes(nNodesPFc, std::vector<double>(i+1));
            std::vector<EMatrix> fcMats(nIPsPFc);
            std::vector<double> fcdV(nIPsPFc, 0.0);
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
            }
            HDGUNabU op(&refEl);
            op.allocate(i+1);
            EVector sol = EVector::Zero(i+1);
            sol[0] = 1.0;
            std::vector<EVector> solution(refEl.getNumNodes(), sol);
            const std::vector< std::vector<double> > * nodes = refEl.getNodes(); 
            for(int k = 0; k < nodes->size(); k++){
              solution[k] *= std::pow(nodes->at(k)[0], n+1);
            }
            CHECK_NOTHROW(op.setSolution(solution));
            std::vector<EVector> trace(nNodesPFc*nFaces, sol);
            for(int iFace = 0; iFace < nFaces; iFace++){
              for(int k = 0; k < nNodesPFc; k++){
                trace[iFace*nNodesPFc + k] = solution[refEl.getFaceNodes()->at(iFace)[k]];
              }
            }
            CHECK_NOTHROW(op.setTrace(trace));
            CHECK_NOTHROW(op.setFromBase(&normals));
            CHECK_NOTHROW(op.assemble(dV, jacs));
            EVector u(refEl.getNumNodes()*(i+1));
            for(int k = 0; k < refEl.getNumNodes(); k++){
              for(int l = 0; l < (i+1); l++){
                u[k*(i+1) + l] = solution[k][l];
              }
            }
            EVector t(nFaces*nNodesPFc*(i+1));
            for(int iFace = 0; iFace < nFaces; iFace++){
              for(int k = 0; k < nNodesPFc; k++){
                for(int l = 0; l < (i+1); l++){
                  t[(iFace*nNodesPFc + k)*(i+1) + l] = trace[iFace*nNodesPFc + k][l];
                }
              }
            }
            double integral = u.transpose()*(op.getMatrix()->block(0, 0, u.size(), u.size()))*u;
            integral += u.transpose()*(op.getMatrix()->block(0, u.size()*(i+2), u.size(), t.size()))*t;
            CHECK(integral/(2.0*(n+1)) == Approx(TestUtils::intx2n(i+1, (3*n+2)/2)).margin(0.2));
            integral = op.getRHS()->segment(0, u.size()).dot(u);
            CHECK(integral/(1.0*(n+1)) == Approx(TestUtils::intx2n(i+1, (3*n+2)/2)).margin(0.2));
          }
        };
      }
    }
  }

};

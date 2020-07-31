#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include <algorithm>
#include "TestUtils.h"
#include "NabUU.h"

using namespace hfox;

TEST_CASE("Testing NabUU operator", "[unit][operator][nonlinear][NabUU]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    NabUU op(&refEl);
    std::vector<EVector> velocity(refEl.getNumNodes(), EVector::Zero(1));
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(op.allocate(1));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(op.setSolution(velocity));
  };

  int maxDim = 3;
  int maxOrd = 5;

  for(int i = 0; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      for(int n = 0; n < (2*(j+1)+1)/3; n++){
        SECTION("Testing NabUU assembly in reference element: (dim = " + std::to_string(i+1) + ", order = " + std::to_string(j+1) + ", n= " + std::to_string(n+1) + ")"){
          if((n+1) % 2 != 0){
            ReferenceElement refEl(i+1, j+1, "simplex");
            NabUU op(&refEl);
            op.allocate(i+1);
            EVector sol = EVector::Zero(i+1);
            sol[0] = 1.0;
            std::vector<EVector> solution(refEl.getNumNodes(), sol);
            const std::vector< std::vector<double> > * nodes = refEl.getNodes(); 
            for(int k = 0; k < nodes->size(); k++){
              solution[k] *= std::pow(nodes->at(k)[0], n+1);
            }
            CHECK_NOTHROW(op.setSolution(solution));
            std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
            CHECK_NOTHROW(op.assemble(*(refEl.getIPWeights()), invJacs));
            EVector u(refEl.getNumNodes()*(i+1));
            for(int k = 0; k < refEl.getNumNodes(); k++){
              for(int l = 0; l < (i+1); l++){
                u[k*(i+1) + l] = solution[k][l];
              }
            }
            double integral = u.transpose()*(*(op.getMatrix()))*u;
            CHECK(integral/(4.0*(n+1)) == Approx(TestUtils::intx2n(i+1, (3*(n+1)-1)/2)));
            integral = op.getRHS()->dot(u);
            CHECK(integral/(2.0*(n+1)) == Approx(TestUtils::intx2n(i+1, (3*(n+1)-1)/2)));
          }
        };
      }
    }
  }

};

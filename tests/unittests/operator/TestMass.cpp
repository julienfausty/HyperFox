#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include <algorithm>
#include "TestUtils.h"
#include "Mass.h"

using namespace hfox;

TEST_CASE("Testing Mass operator", "[unit][operator][Mass]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    Mass massOp(&refEl);
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(massOp.assemble(dummydV, dummyJac));
  };

  std::vector<double> volumes = {2.0, 2.0, 4.0/3.0};

  int maxDim = 3;
  int maxOrd = 5;
  for(int i = 0; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      SECTION("Testing mass in reference element: (dim=" + std::to_string(i+1) + ", order="+std::to_string(j+1)+")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        Mass massOp(&refEl);
        massOp.allocate(1);
        std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
        CHECK_NOTHROW(massOp.assemble(*(refEl.getIPWeights()), invJacs));
        EVector ones = EVector::Constant(refEl.getNumNodes(), 1.0);
        double integral = ones.transpose()*(*(massOp.getMatrix()))*ones;
        CHECK(integral == Approx(volumes[i]));
        std::vector<double> xn(refEl.getNumNodes(), 0.0);
        std::transform(refEl.getNodes()->begin(), refEl.getNodes()->end(), xn.begin(), 
            [j](const std::vector<double> & node){return std::pow(node[0], j+1);});
        EMap<EVector> xnVec(xn.data(), xn.size());
        integral = xnVec.transpose()*(*(massOp.getMatrix()))*xnVec;
        CHECK(integral == Approx(TestUtils::intx2n(i+1, j+1)));
      };

      if(i != 0){
        SECTION("Testing mass in rotated element: (dim=" + std::to_string(i+1) + ", order="+std::to_string(j+1)+")"){
          ReferenceElement refEl(i+1, j+1, "simplex");
          Mass massOp(&refEl);
          massOp.allocate(1);
          std::vector< std::vector<double> > nodes = *(refEl.getNodes());
          EMatrix rot = EMatrix::Identity(i+1, i+1);
          double theta = 1.0;
          EMatrix thetaMat(2, 2);
          thetaMat << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
          rot.block(0,0,2,2) = thetaMat;
          for(int k = 0; k < nodes.size(); k++){
            EMap<EVector> buff(nodes[k].data(), nodes[k].size());
            buff = rot * buff;
          }
          std::vector<EMatrix> jacs = Operator::calcJacobians(nodes, &refEl);
          std::vector<EMatrix> invJacs = Operator::calcInvJacobians(jacs);
          std::vector<double> dV = Operator::calcMeasure(Operator::calcDetJacobians(jacs), &refEl);
          CHECK_NOTHROW(massOp.assemble(dV, invJacs));
          EVector ones = EVector::Constant(refEl.getNumNodes(), 1.0);
          double integral = ones.transpose()*(*(massOp.getMatrix()))*ones;
          CHECK(integral == Approx(volumes[i]));
          std::vector<double> xn(refEl.getNumNodes(), 0.0);
          std::transform(refEl.getNodes()->begin(), refEl.getNodes()->end(), xn.begin(), 
              [j](const std::vector<double> & node){return std::pow(node[0], j+1);});
          EMap<EVector> xnVec(xn.data(), xn.size());
          integral = xnVec.transpose()*(*(massOp.getMatrix()))*xnVec;
          CHECK(integral == Approx(TestUtils::intx2n(i+1, j+1)));
        };
      }
    }
  }

};

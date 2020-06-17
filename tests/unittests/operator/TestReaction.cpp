#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include <algorithm>
#include "TestUtils.h"
#include "Reaction.h"

using namespace hfox;

TEST_CASE("Testing Reaction operator", "[unit][operator][Reaction]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    Reaction op(&refEl);
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(op.allocate(1));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(op.setReactionFunction([](const std::vector<double> & v){return std::pow(v[0], 2);}));
    CHECK_THROWS(op.assemble(dummydV, dummyJac));
  };


  int maxDim = 3;
  int maxOrd = 5;
  for(int i = 0; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      SECTION("Testing reaction in reference element: (dim=" + std::to_string(i+1) + ", order="+std::to_string(j+1)+")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        Reaction op(&refEl);
        op.allocate(1);
        std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
        op.setReactionFunction([j](const std::vector<double> & v){return std::pow(v[0], j+1);});
        CHECK_NOTHROW(op.calcReaction(*(refEl.getNodes())));
        CHECK_NOTHROW(op.assemble(*(refEl.getIPWeights()), invJacs));
        std::vector<double> xn(refEl.getNumNodes(), 0.0);
        std::transform(refEl.getNodes()->begin(), refEl.getNodes()->end(), xn.begin(), 
            [j](const std::vector<double> & node){return std::pow(node[0], j+1);});
        EMap<EVector> xnVec(xn.data(), xn.size()); 
        double integral = EVector::Constant(xn.size(), 1.0).transpose()*(*(op.getMatrix()))*xnVec;
        CHECK(integral == Approx(TestUtils::intx2n(i+1, j+1)));
      };

      if(i != 0){
        SECTION("Testing reaction in rotated element: (dim=" + std::to_string(i+1) + ", order="+std::to_string(j+1)+")"){
          ReferenceElement refEl(i+1, j+1, "simplex");
          Reaction op(&refEl);
          op.allocate(1);
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
          op.setReactionFunction([j](const std::vector<double> & v){return std::pow(v[0], j+1);});
          CHECK_NOTHROW(op.calcReaction(*(refEl.getNodes())));
          CHECK_NOTHROW(op.assemble(dV, invJacs)); 
          std::vector<double> xn(refEl.getNumNodes(), 0.0);
          std::transform(refEl.getNodes()->begin(), refEl.getNodes()->end(), xn.begin(), 
              [j](const std::vector<double> & node){return std::pow(node[0], j+1);});
          EMap<EVector> xnVec(xn.data(), xn.size());
          double integral = EVector::Constant(xn.size(), 1.0).transpose()*(*(op.getMatrix()))*xnVec;
          CHECK(integral == Approx(TestUtils::intx2n(i+1, j+1)));
        };
      }
    }
  }
};

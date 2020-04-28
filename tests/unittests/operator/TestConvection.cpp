#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include <algorithm>
#include "TestUtils.h"
#include "Convection.h"

using namespace hfox;

TEST_CASE("Testing Convection operator", "[unit][operator][Convection]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    Convection convOp(&refEl);
    std::vector<EVector> velocity(refEl.getNumNodes(), EVector::Zero(1));
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(convOp.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(convOp.allocate(1));
    CHECK_THROWS(convOp.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(convOp.setVelocity(velocity));
  };

  int maxDim = 3;
  int maxOrd = 5;
  for(int i = 0; i < maxDim; i++){
    for(int j = 0; j < maxOrd; j++){
      for(int d = 0; d < i+1; d++){
        SECTION("Testing convection assemble in ref element (dim=" + std::to_string(i+1)+ ", ord=" + std::to_string(j+1) + ", e_d=" + std::to_string(d) +")"){
          ReferenceElement refEl(i+1, j+1, "simplex");
          Convection convOp(&refEl);
          convOp.allocate(1);
          EVector vel = EVector::Zero(i+1);
          vel[d] = 1.0;
          std::vector<EVector> velocity(refEl.getNumNodes(), vel);
          CHECK_NOTHROW(convOp.setVelocity(velocity));
          std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
          CHECK_NOTHROW(convOp.assemble(*(refEl.getIPWeights()), invJacs));
          EVector u(refEl.getNumNodes());
          EVector w(refEl.getNumNodes());
          double preFactor = 1.0/(j+1.0);
          for(int k = 0; k < refEl.getNumNodes(); k++){
            u[k] = preFactor*std::pow((refEl.getNodes()->at(k))[d], j+1);
            w[k] = std::pow((refEl.getNodes()->at(k))[d], j);
          }
          double integral = w.transpose()*(*(convOp.getMatrix()))*u;
          CHECK(integral == Approx(TestUtils::intx2n(i+1, j)));
        };

        if(i != 0){
          SECTION("Testing convection assemble in rotated element (dim=" + std::to_string(i+1)+ ", ord=" + std::to_string(j+1) + ", e_d=" + std::to_string(d) +")"){
            ReferenceElement refEl(i+1, j+1, "simplex");
            Convection convOp(&refEl);
            convOp.allocate(1);
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
            EVector vel = EVector::Zero(i+1);
            vel[d] = 1.0;
            vel = rot*vel;
            std::vector<EVector> velocity(refEl.getNumNodes(), vel);
            CHECK_NOTHROW(convOp.setVelocity(velocity));
            CHECK_NOTHROW(convOp.assemble(dV, invJacs));
            EVector u(refEl.getNumNodes());
            EVector w(refEl.getNumNodes());
            double preFactor = 1.0/(j+1.0);
            for(int k = 0; k < refEl.getNumNodes(); k++){
              u[k] = preFactor*std::pow((refEl.getNodes()->at(k))[d], j+1);
              w[k] = std::pow((refEl.getNodes()->at(k))[d], j);
            }
            double integral = w.transpose()*(*(convOp.getMatrix()))*u;
            CHECK(integral == Approx(TestUtils::intx2n(i+1, j)));
          };
        }
      }
    }
  }

};

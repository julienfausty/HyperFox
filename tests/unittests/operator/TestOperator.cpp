#include <catch2/catch.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "DenseEigen.h"
#include "Operator.h"
#include "ReferenceElement.h"
#include "TestUtils.h"

using namespace hfox;

std::vector< std::vector<double> > nonLinearElement(const std::vector< std::vector<double> > & nodes,
    int ord){
  std::vector<double> zeros(nodes[0].size(), 0.0);
  std::vector< std::vector<double> > res(nodes.size(), zeros);
  for(int i = 0; i < nodes.size(); i++){
    for(int j = 0; j < nodes[i].size(); j++){
      res[i][j] = std::pow(nodes[i][j], ord);
    }
  }
  return res;
};

std::vector<EMatrix> nonLinearJacobiansCalc(const std::vector< std::vector<double> > & nodes,
    int ord){
  std::vector<EMatrix> res(nodes.size(), EMatrix::Zero(nodes[0].size(), nodes[0].size()));
  for(int i = 0; i < nodes.size(); i++){
    for(int j = 0; j < nodes[i].size(); j++){
      res[i](j, j) = ord*std::pow(nodes[i][j], ord-1);
    }
  }
  return res;
};

TEST_CASE("Testing static methods in Operator", "[unit][operator][Operator]"){
  EMatrix Jac(2, 2);
  Jac << 3, 2,
      1, 4;
  std::vector<double> offset(2, 0.0);
  ReferenceElement * refEl = new ReferenceElement(2, 3, "simplex");
  std::vector< std::vector<double> > linEl = TestUtils::linElement(*(refEl->getNodes()), Jac, offset);
  SECTION("Test linear element Jacobian computation (2D)"){
    std::vector<EMatrix> jacobians = Operator::calcJacobians(linEl, refEl);
    std::vector<EMatrix> invJacobians = Operator::calcInvJacobians(jacobians);
    std::vector<double> detJacobians = Operator::calcDetJacobians(jacobians);
    EMatrix invJT = Jac.transpose().inverse();
    double detJT = Jac.transpose().determinant();
    for(int i = 0; i < jacobians.size(); i++){
      for(int j = 0; j < jacobians[i].rows(); j++){
        for(int k = 0; k < jacobians[i].cols(); k++){
          CHECK(jacobians[i](j, k) == Approx(Jac(k, j)).margin(1e-12));
          CHECK(invJacobians[i](j, k) == Approx(invJT(j, k)).margin(1e-12));
        }
      }
      CHECK(detJacobians[i] == Approx(detJT).margin(1e-12));
    }
  };

  Jac(0,0) += 2; Jac(0,1) += -1;
  Jac(1, 0) += 1; Jac(1,1) += 3;
  offset[0] = 1.0; offset[1] = 2.0;
  linEl = TestUtils::linElement(*(refEl->getNodes()), Jac, offset);
  SECTION("Test linear element Jacobian computation with offset (2D)"){
    std::vector<EMatrix> jacobians = Operator::calcJacobians(linEl, refEl);
    std::vector<EMatrix> invJacobians = Operator::calcInvJacobians(jacobians);
    std::vector<double> detJacobians = Operator::calcDetJacobians(jacobians);
    EMatrix invJT = Jac.transpose().inverse();
    double detJT = Jac.transpose().determinant();
    for(int i = 0; i < jacobians.size(); i++){
      for(int j = 0; j < jacobians[i].rows(); j++){
        for(int k = 0; k < jacobians[i].cols(); k++){
          CHECK(jacobians[i](j, k) == Approx(Jac(k, j)).margin(1e-12));
          CHECK(invJacobians[i](j, k) == Approx(invJT(j, k)).margin(1e-12));
        }
      }
      CHECK(detJacobians[i] == Approx(detJT).margin(1e-12));
    }
  };

  std::vector< std::vector<double> > nonLinEl = nonLinearElement(*(refEl->getNodes()), 2);
  const std::vector< std::vector<double> > refNodes = *(refEl->getNodes());
  std::vector<EMatrix> anaJacs = nonLinearJacobiansCalc(*(refEl->getIPCoords()), 2);
  SECTION("Test non linear element Jacobian computation (2D)"){
    std::vector<EMatrix> jacobians = Operator::calcJacobians(nonLinEl, refEl);
    std::vector<EMatrix> invJacobians = Operator::calcInvJacobians(jacobians);
    std::vector<double> detJacobians = Operator::calcDetJacobians(jacobians);
    CHECK(jacobians.size() == anaJacs.size());
    for(int i = 0; i < jacobians.size(); i++){
      EMatrix invJT = anaJacs[i].transpose().inverse();
      double detJT = anaJacs[i].transpose().determinant();
      for(int j = 0; j < jacobians[i].rows(); j++){
        for(int k = 0; k < jacobians[i].cols(); k++){
          CHECK(jacobians[i](j, k) == Approx(anaJacs[i](k, j)).margin(1e-12));
          CHECK(invJacobians[i](j, k) == Approx(invJT(j, k)).margin(1e-12));
        }
      }
      CHECK(detJacobians[i] == Approx(detJT).margin(1e-12));
    }
  };

  delete refEl;
  refEl = new ReferenceElement(3, 4, "simplex");
  Jac.resize(3, 3);
  Jac(0, 0) = 1;
  Jac(2, 2) = 5; 
  offset.resize(3, 0.0);
  linEl = TestUtils::linElement(*(refEl->getNodes()), Jac, offset);
  SECTION("Test linear element Jacobian computation (3D)"){
    std::vector<EMatrix> jacobians = Operator::calcJacobians(linEl, refEl);
    std::vector<EMatrix> invJacobians = Operator::calcInvJacobians(jacobians);
    std::vector<double> detJacobians = Operator::calcDetJacobians(jacobians);
    EMatrix invJT = Jac.transpose().inverse();
    double detJT = Jac.transpose().determinant();
    for(int i = 0; i < jacobians.size(); i++){
      for(int j = 0; j < jacobians[i].rows(); j++){
        for(int k = 0; k < jacobians[i].cols(); k++){
          CHECK(jacobians[i](j, k) == Approx(Jac(k, j)).margin(1e-12));
          CHECK(invJacobians[i](j, k) == Approx(invJT(j, k)).margin(1e-12));
        }
      }
      CHECK(detJacobians[i] == Approx(detJT).margin(1e-12));
    }
  };

  Jac(0,0) += 2; Jac(0,1) += -1;
  Jac(1, 0) += 1; Jac(1,1) += 3;
  offset[0] = 1.0; offset[1] = 2.0;
  linEl = TestUtils::linElement(*(refEl->getNodes()), Jac, offset);
  SECTION("Test linear element Jacobian computation with offset (3D)"){
    std::vector<EMatrix> jacobians = Operator::calcJacobians(linEl, refEl);
    std::vector<EMatrix> invJacobians = Operator::calcInvJacobians(jacobians);
    std::vector<double> detJacobians = Operator::calcDetJacobians(jacobians);
    EMatrix invJT = Jac.transpose().inverse();
    double detJT = Jac.transpose().determinant();
    for(int i = 0; i < jacobians.size(); i++){
      for(int j = 0; j < jacobians[i].rows(); j++){
        for(int k = 0; k < jacobians[i].cols(); k++){
          CHECK(jacobians[i](j, k) == Approx(Jac(k, j)).margin(1e-12));
          CHECK(invJacobians[i](j, k) == Approx(invJT(j, k)).margin(1e-12));
        }
      }
      CHECK(detJacobians[i] == Approx(detJT).margin(1e-12));
    }
  };

  nonLinEl = nonLinearElement(*(refEl->getNodes()), 3);
  anaJacs = nonLinearJacobiansCalc(*(refEl->getIPCoords()), 3);
  SECTION("Test non linear element Jacobian computation (3D)"){
    std::vector<EMatrix> jacobians = Operator::calcJacobians(nonLinEl, refEl);
    std::vector<EMatrix> invJacobians = Operator::calcInvJacobians(jacobians);
    std::vector<double> detJacobians = Operator::calcDetJacobians(jacobians);
    CHECK(jacobians.size() == anaJacs.size());
    for(int i = 0; i < jacobians.size(); i++){
      EMatrix invJT = anaJacs[i].transpose().inverse();
      double detJT = anaJacs[i].transpose().determinant();
      for(int j = 0; j < jacobians[i].rows(); j++){
        for(int k = 0; k < jacobians[i].cols(); k++){
          CHECK(jacobians[i](j, k) == Approx(anaJacs[i](k, j)).margin(1e-12));
          CHECK(invJacobians[i](j, k) == Approx(invJT(j, k)).margin(1e-12));
        }
      }
      CHECK(detJacobians[i] == Approx(detJT).margin(1e-12));
    }
  };

  delete refEl;

};

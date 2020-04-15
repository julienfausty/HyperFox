#include <catch2/catch.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include "DenseEigen.h"
#include "Operator.h"
#include "Diffusion.h"
#include "ReferenceElement.h"
#include "TestUtils.h"

using namespace hfox;

double monomials(const std::vector<double> & point, int ord){
  std::vector<double> monos(point.size());
  std::transform(point.begin(), point.end(), monos.begin(), [ord](double x){return std::pow(x, ord);});
  return std::accumulate(monos.begin(), monos.end(), 0.0);
}

TEST_CASE("Testing the Diffusion operator", "[unit][operator][Diffusion]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    Diffusion diffOp(&refEl);
    std::vector<double> dummyDet(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(diffOp.assemble(dummyDet, dummyJac));
  };

  int maxDim = 2, maxOrder = 5;
  ReferenceElement * refEl;
  std::vector<double> integrals(maxDim*maxOrder);
  for(int j = 0; j < maxOrder; j++){
    integrals[j] = 2.0*std::pow(j+1, 2.0)/(2.0*(j+1)-1.0);
    integrals[maxOrder + j] = 2.0*integrals[j];
  }
  //dim loop
  for(int i = 0; i < maxDim; i++){
    //order loop
    for(int j = 0; j < maxOrder; j++){
      SECTION("Test Diffusion operator on monomials in ref element (dim=" + std::to_string(i+1) + ", order" + std::to_string(j+1) + ")"){
        refEl = new ReferenceElement(i+1, j+1, "simplex");
        Diffusion diffOp(refEl);
        diffOp.allocate(1);
        std::vector<EMatrix> jacobians = Operator::calcJacobians(*(refEl->getNodes()), refEl);
        std::vector<double> detJacs(jacobians.size());
        std::vector<EMatrix> invJacs(jacobians.size());
        std::transform(jacobians.begin(), jacobians.end(), detJacs.begin(), [](EMatrix jac){return jac.determinant();});
        std::transform(jacobians.begin(), jacobians.end(), invJacs.begin(), [](EMatrix jac){return jac.inverse();});
        CHECK_NOTHROW(diffOp.assemble(detJacs, invJacs));
        std::vector< double > monos(refEl->getNumNodes());
        std::transform(refEl->getNodes()->begin(), refEl->getNodes()->end(), monos.begin(), [j](const std::vector<double> & point){return monomials(point, j+1);});
        Eigen::Map<EVector> vecMonos(monos.data(), monos.size());
        double integral = vecMonos.transpose() * (*(diffOp.getMatrix())) * vecMonos;
        CHECK(integral == Approx(integrals[i*maxOrder + j]).margin(1e-12));
        delete refEl;
      };
    }
  }


  //dim loop
  for(int i = 1; i < maxDim; i++){
    //order loop
    for(int j = 0; j < maxOrder; j++){
      SECTION("Test Diffusion operator on monomials in rotated element (dim=" + std::to_string(i+1) + ", order" + std::to_string(j+1) + ")"){
        refEl = new ReferenceElement(i+1, j+1, "simplex");
        Diffusion diffOp(refEl);
        diffOp.allocate(1);
        EMatrix elJac = EMatrix::Identity(i+1, i+1); 
        double theta = 1.0;//rotation angle
        EMatrix rot(2,2); rot << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
        elJac.block<2,2>(0,0) = rot;
        double elDet = elJac.determinant();
        std::vector<double> offset(i+1, 2.0);
        std::vector< std::vector<double> > linElement = TestUtils::linElement(*(refEl->getNodes()), elJac, offset);
        std::vector<EMatrix> jacobians = Operator::calcJacobians(linElement, refEl);
        std::vector<double> detJacs(jacobians.size());
        std::vector<EMatrix> invJacs(jacobians.size());
        std::transform(jacobians.begin(), jacobians.end(), detJacs.begin(), [](EMatrix jac){return jac.determinant();});
        std::transform(jacobians.begin(), jacobians.end(), invJacs.begin(), [](EMatrix jac){return jac.inverse();});
        CHECK_NOTHROW(diffOp.assemble(detJacs, invJacs));
        std::vector< double > monos(refEl->getNumNodes());
        std::transform(refEl->getNodes()->begin(), refEl->getNodes()->end(), monos.begin(), [j](const std::vector<double> & point){return monomials(point, j+1);});
        Eigen::Map<EVector> vecMonos(monos.data(), monos.size());
        double integral = vecMonos.transpose() * (*(diffOp.getMatrix())) * vecMonos;
        CHECK(integral == Approx(integrals[i*maxOrder + j]).margin(1e-12));
        delete refEl;
      };
    }
  }

  //dim loop
  for(int i = 1; i < maxDim; i++){
    //order loop
    for(int j = 0; j < maxOrder; j++){
      SECTION("Test Diffusion operator on monomials in rotated element with multiple allocation (dim=" + std::to_string(i+1) + ", order" + std::to_string(j+1) + ")"){
        refEl = new ReferenceElement(i+1, j+1, "simplex");
        Diffusion diffOp(refEl);
        diffOp.allocate(2);
        EMatrix elJac = EMatrix::Identity(i+1, i+1); 
        double theta = 1.0;//rotation angle
        EMatrix rot(2,2); rot << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
        elJac.block<2,2>(0,0) = rot;
        double elDet = elJac.determinant();
        std::vector<double> offset(i+1, 2.0);
        std::vector< std::vector<double> > linElement = TestUtils::linElement(*(refEl->getNodes()), elJac, offset);
        std::vector<EMatrix> jacobians = Operator::calcJacobians(linElement, refEl);
        std::vector<double> detJacs(jacobians.size());
        std::vector<EMatrix> invJacs(jacobians.size());
        std::transform(jacobians.begin(), jacobians.end(), detJacs.begin(), [](EMatrix jac){return jac.determinant();});
        std::transform(jacobians.begin(), jacobians.end(), invJacs.begin(), [](EMatrix jac){return jac.inverse();});
        CHECK_NOTHROW(diffOp.assemble(detJacs, invJacs));
        int nNodes = refEl->getNumNodes();
        std::vector< double > monos(nNodes);
        std::transform(refEl->getNodes()->begin(), refEl->getNodes()->end(), monos.begin(), [j](const std::vector<double> & point){return monomials(point, j+1);});
        std::vector<double> doubleMonos(2*nNodes);
        for(int k = 0; k < nNodes; k++){
          doubleMonos[2*k] = monos[k];
          doubleMonos[2*k + 1] = monos[k];
        }
        Eigen::Map<EVector> vecMonos(doubleMonos.data(), doubleMonos.size());
        double integral = vecMonos.transpose() * (*(diffOp.getMatrix())) * vecMonos;
        CHECK(integral == Approx(2.0*integrals[i*maxOrder + j]).margin(1e-12));
        delete refEl;
      };
    }
  }

  //dim loop
  for(int i = 1; i < maxDim; i++){
    //order loop
    for(int j = 0; j < maxOrder; j++){
      SECTION("Test Diffusion operator on monomials in rotated element with scalar diffusion (dim=" + std::to_string(i+1) + ", order" + std::to_string(j+1) + ")"){
        refEl = new ReferenceElement(i+1, j+1, "simplex");
        Diffusion diffOp(refEl);
        diffOp.allocate(1);
        EMatrix elJac = EMatrix::Identity(i+1, i+1); 
        double theta = 1.5;//rotation angle
        EMatrix rot(2,2); rot << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
        elJac.block<2,2>(0,0) = rot;
        double elDet = elJac.determinant();
        std::vector<double> offset(i+1, 2.0);
        std::vector< std::vector<double> > linElement = TestUtils::linElement(*(refEl->getNodes()), elJac, offset);
        std::vector<EMatrix> jacobians = Operator::calcJacobians(linElement, refEl);
        std::vector<double> detJacs(jacobians.size());
        std::vector<EMatrix> invJacs(jacobians.size());
        std::transform(jacobians.begin(), jacobians.end(), detJacs.begin(), [](EMatrix jac){return jac.determinant();});
        std::transform(jacobians.begin(), jacobians.end(), invJacs.begin(), [](EMatrix jac){return jac.inverse();});
        std::vector<EMatrix> diffCoeff(refEl->getNumNodes(), (1.0/(i+1.0))*EMatrix::Identity(1, 1));
        CHECK_NOTHROW(diffOp.setDiffTensor(diffCoeff));
        CHECK_NOTHROW(diffOp.assemble(detJacs, invJacs));
        std::vector< double > monos(refEl->getNumNodes());
        std::transform(refEl->getNodes()->begin(), refEl->getNodes()->end(), monos.begin(), [j](const std::vector<double> & point){return monomials(point, j+1);});
        Eigen::Map<EVector> vecMonos(monos.data(), monos.size());
        double integral = vecMonos.transpose() * (*(diffOp.getMatrix())) * vecMonos;
        CHECK(integral == Approx((1.0/(i+1.0))*integrals[i*maxOrder + j]).margin(1e-12));
        delete refEl;
      };
    }
  }

  //dim loop
  for(int i = 1; i < maxDim; i++){
    //order loop
    for(int j = 0; j < maxOrder; j++){
      SECTION("Test Diffusion operator on monomials in rotated element with tensorial diffusion (dim=" + std::to_string(i+1) + ", order" + std::to_string(j+1) + ")"){
        refEl = new ReferenceElement(i+1, j+1, "simplex");
        Diffusion diffOp(refEl);
        diffOp.allocate(1);
        EMatrix elJac = EMatrix::Identity(i+1, i+1); 
        double theta = 1.5;//rotation angle
        EMatrix rot(2,2); rot << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
        elJac.block<2,2>(0,0) = rot;
        double elDet = elJac.determinant();
        std::vector<double> offset(i+1, 2.0);
        std::vector< std::vector<double> > linElement = TestUtils::linElement(*(refEl->getNodes()), elJac, offset);
        std::vector<EMatrix> jacobians = Operator::calcJacobians(linElement, refEl);
        std::vector<double> detJacs(jacobians.size());
        std::vector<EMatrix> invJacs(jacobians.size());
        std::transform(jacobians.begin(), jacobians.end(), detJacs.begin(), [](EMatrix jac){return jac.determinant();});
        std::transform(jacobians.begin(), jacobians.end(), invJacs.begin(), [](EMatrix jac){return jac.inverse();});
        std::vector<EMatrix> diffCoeff(refEl->getNumNodes(), (1.0/(i+1.0))*EMatrix::Identity(i+1, i+1));
        CHECK_NOTHROW(diffOp.setDiffTensor(diffCoeff));
        CHECK_NOTHROW(diffOp.assemble(detJacs, invJacs));
        std::vector< double > monos(refEl->getNumNodes());
        std::transform(refEl->getNodes()->begin(), refEl->getNodes()->end(), monos.begin(), [j](const std::vector<double> & point){return monomials(point, j+1);});
        Eigen::Map<EVector> vecMonos(monos.data(), monos.size());
        double integral = vecMonos.transpose() * (*(diffOp.getMatrix())) * vecMonos;
        CHECK(integral == Approx((1.0/(i+1.0))*integrals[i*maxOrder + j]).margin(1e-12));
        delete refEl;
      };
    }
  }

};

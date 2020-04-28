#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include <algorithm>
#include "TestUtils.h"
#include "HDGBase.h"

using namespace hfox;

std::vector< std::vector<double> > getRefNormals(int dim){
  std::vector< std::vector<double> > normals;
  if(dim == 2){
    normals.push_back(std::vector<double>({0, -1}));
    normals.push_back(std::vector<double>({1.0/std::sqrt(2.0), 1.0/std::sqrt(2.0)}));
    normals.push_back(std::vector<double>({-1, 0}));
  } else if(dim == 3){
    normals.push_back(std::vector<double>({0, -1, 0}));
    normals.push_back(std::vector<double>({1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)}));
    normals.push_back(std::vector<double>({-1, 0, 0}));
    normals.push_back(std::vector<double>({0, 0, -1}));
  }
  return normals;
};//getRefNormals

TEST_CASE("Testing HDGBase operator", "[unit][operator][HDGBase]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(1, 1, "simplex");
    CHECK_NOTHROW(HDGBase(&refEl));
    HDGBase baseOp(&refEl);
    std::vector<double> tau(refEl.getFaceElement()->getNumNodes() * refEl.getNumFaces(), 0.0);
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(baseOp.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(baseOp.allocate(1));
    CHECK_THROWS(baseOp.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(baseOp.setTau(tau));
  };

  int maxDim = 3;
  int maxOrder = 5;

  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrder; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      HDGBase baseOp(&refEl);
      baseOp.allocate(1);
      std::vector<double> taus(refEl.getFaceElement()->getNumNodes() * refEl.getNumFaces(), 1.0);
      baseOp.setTau(taus);
      std::vector<EMatrix> invJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
      baseOp.assemble(*(refEl.getIPWeights()), invJacs);
      SECTION("Testing bulk operator parts of base in reference element (dim=" + std::to_string(i+1)+ ", ord=" + std::to_string(j+1) + ")"){
        int nNodes = refEl.getNumNodes();
        Mass bulkMassOp(&refEl);
        //Sqq
        bulkMassOp.allocate(i+1);
        bulkMassOp.assemble(*(refEl.getIPWeights()), invJacs);
        CHECK(baseOp.getMatrix()->block(nNodes, nNodes, (i+1)*nNodes, (i+1)*nNodes) == *(bulkMassOp.getMatrix()));
        //Squ
        Convection convOp(&refEl);
        convOp.allocate(1);
        EMatrix Squ(nNodes*(i+1), nNodes);
        for(int k = 0; k < i+1; k++){
          EVector unitBase = EVector::Zero(i+1);
          unitBase[k] = 1.0;
          std::vector<EVector> velocity(nNodes, unitBase);
          convOp.setVelocity(velocity);
          convOp.assemble(*(refEl.getIPWeights()), invJacs);
          for(int l = 0 ; l < nNodes; l++){
            Squ.row(l*(i+1) + k) = (convOp.getMatrix()->col(l)).transpose();
          }
        }
        CHECK(baseOp.getMatrix()->block(nNodes, 0, nNodes*(i+2), nNodes) == Squ);
      };
      SECTION("Testing surface operator parts of base in reference element (dim=" + std::to_string(i+1)+ ", ord=" + std::to_string(j+1) + ")"){
        const ReferenceElement * fEl = refEl.getFaceElement();
        Mass faceMassOp(fEl);
        faceMassOp.allocate(1);
        int nNodes = refEl.getNumNodes();
        int nNodesFace = fEl->getNumNodes();
        std::vector<EMatrix> faceOps(refEl.getNumFaces(), EMatrix::Zero(nNodesFace, nNodesFace));
        std::vector<EMatrix> jacs(fEl->getNumIPs(), EMatrix::Zero(i, i+1));
        std::vector<EMatrix> invJacs(fEl->getNumIPs(), EMatrix::Zero(i+1, i));
        std::vector<double> dV(fEl->getNumIPs(), 0.0);
        const std::vector< std::vector<int> > * faceNodeMap = refEl.getFaceNodes();
        std::vector< std::vector<double> > faceNodes(fEl->getNumNodes(), std::vector<double>(i+1, 0.0));
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          for(int k = 0; k < fEl->getNumNodes(); k++){
            faceNodes[k] = refEl.getNodes()->at((faceNodeMap->at(iface))[k]);
          }
          jacs = Operator::calcJacobians(faceNodes, fEl);
          invJacs = Operator::calcInvJacobians(jacs);
          dV = Operator::calcMeasure(Operator::calcDetJacobians(jacs), fEl);
          faceMassOp.assemble(dV, invJacs);
          faceOps[iface] = *(faceMassOp.getMatrix());
          int startDiag = nNodes*(i+2) + nNodesFace*iface;
          //Sll
          CHECK(baseOp.getMatrix()->block(startDiag, startDiag, nNodesFace, nNodesFace) == -faceOps[iface]);
        }
        //Slu
        EMatrix Slu = EMatrix::Zero(nNodesFace*refEl.getNumFaces(), nNodes);
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          for(int k = 0; k < nNodesFace; k++){
            Slu.col((faceNodeMap->at(iface))[k]) += (faceOps[iface]).col(k);
          }
        }
        CHECK(baseOp.getMatrix()->block(nNodes*(i+2), 0, nNodesFace*refEl.getNumFaces(), nNodes) == Slu);
        //Sul
        CHECK(baseOp.getMatrix()->block(0, nNodes*(i+2), nNodes, nNodesFace*refEl.getNumFaces()) == Slu.transpose());
        //Suu
        EMatrix Suu = EMatrix::Zero(nNodes, nNodes);
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          for(int k = 0; k < nNodesFace; k++){
            Suu.row((faceNodeMap->at(iface))[k]) += Slu.row(iface*nNodesFace + k);
          }
        }
        CHECK(baseOp.getMatrix()->block(0, 0, nNodes, nNodes) == -Suu);
        //Slq
        std::vector< std::vector<double> > normals = getRefNormals(i+1);
        EMatrix Slq = EMatrix::Zero(nNodesFace*refEl.getNumFaces(), nNodes*(i+1));
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          for(int k = 0; k < nNodesFace; k++){
            for(int d = 0; d < (i+1); d++){
              Slq.col(faceNodeMap->at(iface)[k]*(i+1)+d) = faceOps[iface].col(k)*normals[iface][d];
            }
          }
        }
        CHECK(baseOp.getMatrix()->block(nNodes*(i+2), nNodes, Slq.rows(), Slq.cols()) == Slq);
        //Sql
        CHECK(baseOp.getMatrix()->block(nNodes, nNodes*(i+2), Slq.cols(), Slq.rows()) == -Slq.transpose());
        //Suq
        CHECK(baseOp.getMatrix()->block(0, nNodes, nNodes, nNodes*(i+1)) == EMatrix::Zero(nNodes, nNodes*(i+1)));
      };
    }
  }
};

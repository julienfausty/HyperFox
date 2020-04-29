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
    std::vector<double> dummydV(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(baseOp.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(baseOp.allocate(1));
    CHECK_THROWS(baseOp.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(baseOp.setTau(tau));
    CHECK_THROWS(baseOp.assemble(dummydV, dummyJac));
  };

  int maxDim = 2;
  int maxOrder = 1;
  for(int i = 1; i < maxDim; i++){
    for(int j = 0; j < maxOrder; j++){
      ReferenceElement refEl(i+1, j+1, "simplex");
      HDGBase baseOp(&refEl);
      baseOp.allocate(1);
      std::vector<double> taus(refEl.getFaceElement()->getNumNodes() * refEl.getNumFaces(), 1.0);
      baseOp.setTau(taus);
      std::vector<EMatrix> elJacs(refEl.getNumIPs(), EMatrix::Identity(i+1, i+1));
      std::vector<EMatrix> jacs(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs());
      std::vector<EMatrix> invJacs(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs());
      std::vector<double> dV(refEl.getNumIPs()+refEl.getNumFaces()*refEl.getFaceElement()->getNumIPs(), 0.0);
      std::copy(elJacs.begin(), elJacs.end(), jacs.begin());
      std::copy(elJacs.begin(), elJacs.end(), invJacs.begin());
      std::copy(refEl.getIPWeights()->begin(), refEl.getIPWeights()->end(), dV.begin());
      std::vector<EMatrix> faceJacs(refEl.getFaceElement()->getNumIPs(), EMatrix::Zero(i, i+1));
      std::vector<double> facedV(refEl.getFaceElement()->getNumIPs(), 0.0);
      std::vector< std::vector<double> > faceNodes(refEl.getFaceElement()->getNumNodes(), std::vector<double>(i+1, 0.0));
      for(int k = 0; k < refEl.getNumFaces(); k++){
        for(int l = 0; l < refEl.getFaceElement()->getNumNodes(); l++){
          faceNodes[l] = refEl.getNodes()->at(refEl.getFaceNodes()->at(k)[l]);
        }
        int offset = refEl.getNumIPs() + k*(refEl.getFaceElement()->getNumIPs());
        faceJacs = Operator::calcJacobians(faceNodes, refEl.getFaceElement());
        std::copy(faceJacs.begin(), faceJacs.end(), jacs.begin() + offset);
        facedV = Operator::calcMeasure(Operator::calcDetJacobians(faceJacs), refEl.getFaceElement());
        std::copy(facedV.begin(), facedV.end(), dV.begin() + offset);
        faceJacs = Operator::calcInvJacobians(faceJacs);
        std::copy(faceJacs.begin(), faceJacs.end(), invJacs.begin() + offset);
      }
      CHECK_NOTHROW(baseOp.calcNormals(*(refEl.getNodes()), jacs));
      baseOp.assemble(dV, invJacs);
      SECTION("Testing bulk operator parts of base in reference element (dim=" + std::to_string(i+1)+ ", ord=" + std::to_string(j+1) + ")"){
        int nNodes = refEl.getNumNodes();
        Mass bulkMassOp(&refEl);
        //Sqq
        bulkMassOp.allocate(i+1);
        bulkMassOp.assemble(*(refEl.getIPWeights()), invJacs);
        CHECK(std::abs((baseOp.getMatrix()->block(nNodes, nNodes, (i+1)*nNodes, (i+1)*nNodes) - *(bulkMassOp.getMatrix())).sum()) < 1e-12);
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
        CHECK(std::abs((baseOp.getMatrix()->block(nNodes, 0, nNodes*(i+1), nNodes) - Squ).sum()) < 1e-12);
      };
      SECTION("Testing surface operator parts of base in reference element (dim=" + std::to_string(i+1)+ ", ord=" + std::to_string(j+1) + ")"){
        const ReferenceElement * fEl = refEl.getFaceElement();
        Mass faceMassOp(fEl);
        faceMassOp.allocate(1);
        int nNodes = refEl.getNumNodes();
        int nNodesFace = fEl->getNumNodes();
        std::vector<EMatrix> faceOps(refEl.getNumFaces(), EMatrix::Zero(nNodesFace, nNodesFace));
        const std::vector< std::vector<int> > * faceNodeMap = refEl.getFaceNodes();
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          int offset = refEl.getNumIPs() + iface*(refEl.getFaceElement()->getNumIPs());
          std::copy(invJacs.begin() + offset, invJacs.begin() + offset + refEl.getFaceElement()->getNumIPs(), faceJacs.begin());
          std::copy(dV.begin() + offset, dV.begin() + offset + refEl.getFaceElement()->getNumIPs(), facedV.begin());
          faceMassOp.assemble(facedV, faceJacs);
          faceOps[iface] = *(faceMassOp.getMatrix());
          int startDiag = nNodes*(i+2) + nNodesFace*iface;
          //Sll 
          CHECK(std::abs((baseOp.getMatrix()->block(startDiag, startDiag, nNodesFace, nNodesFace) + faceOps[iface]).sum()) < 1e-12);
        }
        //Slu
        EMatrix Slu = EMatrix::Zero(nNodesFace*refEl.getNumFaces(), nNodes);
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          for(int k = 0; k < nNodesFace; k++){
            Slu.block(iface*nNodesFace, (faceNodeMap->at(iface))[k], nNodesFace, 1) += (faceOps[iface]).col(k);
          }
        }
        CHECK(std::abs((baseOp.getMatrix()->block(nNodes*(i+2), 0, nNodesFace*refEl.getNumFaces(), nNodes) - Slu).sum()) < 1e-12);
        //Sul
        CHECK(std::abs((baseOp.getMatrix()->block(0, nNodes*(i+2), nNodes, nNodesFace*refEl.getNumFaces()) - Slu.transpose()).sum()) < 1e-12);
        //Suu
        EMatrix Suu = EMatrix::Zero(nNodes, nNodes);
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          for(int k = 0; k < nNodesFace; k++){
            Suu.row((faceNodeMap->at(iface))[k]) += Slu.row(iface*nNodesFace + k);
          }
        }
        CHECK(std::abs((baseOp.getMatrix()->block(0, 0, nNodes, nNodes) + Suu).sum()) < 1e-12);
        //Slq
        std::vector< std::vector<double> > normals = getRefNormals(i+1);
        EMatrix Slq = EMatrix::Zero(nNodesFace*refEl.getNumFaces(), nNodes*(i+1));
        for(int iface = 0; iface < refEl.getNumFaces(); iface++){
          for(int k = 0; k < nNodesFace; k++){
            for(int d = 0; d < (i+1); d++){
              Slq.block(iface*nNodesFace, faceNodeMap->at(iface)[k]*(i+1)+d, nNodesFace, 1) = faceOps[iface].col(k)*normals[iface][d];
            }
          }
        }
        CHECK(std::abs((baseOp.getMatrix()->block(nNodes*(i+2), nNodes, Slq.rows(), Slq.cols()) - Slq).sum()) < 1e-12);
        //Sql
        CHECK(std::abs((baseOp.getMatrix()->block(nNodes, nNodes*(i+2), Slq.cols(), Slq.rows()) + Slq.transpose()).sum()) < 1e-12);
        //Suq
        CHECK(std::abs((baseOp.getMatrix()->block(0, nNodes, nNodes, nNodes*(i+1)) - EMatrix::Zero(nNodes, nNodes*(i+1))).sum()) < 1e-12);
      };
    }
  }
};

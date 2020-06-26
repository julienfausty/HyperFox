#ifndef TESTUTILS_H
#define TESTUTILS_H

#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <chrono>
#include <functional>
#include <cmath>
#define _USE_MATH_DEFINES
#include "DenseEigen.h"
#include "Field.h"
#include "Mesh.h"
#include "Operator.h"


namespace hfox{

class TestUtils{

  public:

  static std::string getRessourcePath(){
    std::string cmake_src = "/home/jfausty/workspace/Codes/HyperFox";
    return std::string(cmake_src + "/ressources");
  };//getRessourcePath

  static bool checkFilesEqual(const std::string& lFilePath, const std::string& rFilePath){
    std::ifstream lFile(lFilePath.c_str(), std::ifstream::in | std::ifstream::binary);
    std::ifstream rFile(rFilePath.c_str(), std::ifstream::in | std::ifstream::binary);

    if(!lFile.is_open() || !rFile.is_open())
    {
      return false;
    }

    int buff_size = 1024;

    char *lBuffer = new char[buff_size]();
    char *rBuffer = new char[buff_size]();

    do {
      lFile.read(lBuffer, buff_size);
      rFile.read(rBuffer, buff_size);

      if (std::memcmp(lBuffer, rBuffer, buff_size) != 0)
      {
        delete[] lBuffer;
        delete[] rBuffer;
        return false;
      }
    } while (lFile.good() || rFile.good());

    delete[] lBuffer;
    delete[] rBuffer;
    lFile.close();
    rFile.close();
    return true;
  };

  template <class T>
    static std::vector<T> unpack(std::vector< std::vector<T> > v){
      std::vector<T> res(v.size()*v[0].size());
      int index = 0;
      for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
          res[index] = v[i][j];
          index += 1;
        }
      }
      return res;
    };
  
  static std::vector< std::vector<double> > linElement(const std::vector< std::vector<double> > & nodes, 
      const EMatrix & Jac, const std::vector<double> & offset){
    std::vector< std::vector<double> > res(nodes.size(), std::vector<double>(Jac.rows(), 0.0));
    Eigen::Map<const EVector> shift(offset.data(), offset.size());
    for(int i = 0; i < nodes.size(); i++){
      Eigen::Map<const EVector> refNode(nodes[i].data(), nodes[i].size()); 
      Eigen::Map<EVector> elNode(res[i].data(), res[i].size());
      elNode = Jac*refNode + shift;
    }
    return res;
  };

  static double intx2n(int dim, int ord){
    double res = 0;
    if(dim < 3){
      res = 2.0/(2.0*ord + 1);
    } else if(dim == 3){
      res = 1.0/(2.0*ord+3.0) + 1.0/(2.0*ord + 1.0);
    }
    return res;
  };
  
  static double l2ProjectionNodeField(Field * field1, Field * field2, Mesh * mesh){
    double integral = 0;
    const ReferenceElement * refEl = mesh->getReferenceElement();
    std::vector<int> cell(refEl->getNumNodes(), 0);
    std::vector< std::vector<double> > nodes(cell.size(), 
        std::vector<double>(mesh->getNodeSpaceDimension(), 0.0));
    const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
    const std::vector<double> * ipWeights = refEl->getIPWeights();
    std::vector<double> detJacs(refEl->getNumIPs(), 0.0);
    for(int i = 0; i < mesh->getNumberCells(); i++){
      mesh->getCell(i, &cell);
      mesh->getSlicePoints(cell, &nodes);
      detJacs = Operator::calcDetJacobians(Operator::calcJacobians(nodes, refEl));
      for(int j = 0; j < refEl->getNumIPs(); j++){
        double dV = detJacs[j]*ipWeights->at(j);
        for(int k = 0; k < cell.size(); k++){
          for(int l = 0; l < cell.size(); l++){
            integral += ipShapes->at(j)[l]*(ipShapes->at(j)[k])*(field1->getValues()->at(cell[k]))*(field2->getValues()->at(cell[l]))*dV;
          }
        }
      }
    }
    return integral;
  };

  static double l2ProjectionCellField(Field * field1, Field * field2, Mesh * mesh){
    double integral = 0;
    const ReferenceElement * refEl = mesh->getReferenceElement();
    std::vector<int> cell(refEl->getNumNodes(), 0);
    std::vector< std::vector<double> > nodes(cell.size(), 
        std::vector<double>(mesh->getNodeSpaceDimension(), 0.0));
    const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
    const std::vector<double> * ipWeights = refEl->getIPWeights();
    std::vector<double> detJacs(refEl->getNumIPs(), 0.0);
    for(int i = 0; i < mesh->getNumberCells(); i++){
      mesh->getCell(i, &cell);
      mesh->getSlicePoints(cell, &nodes);
      detJacs = Operator::calcDetJacobians(Operator::calcJacobians(nodes, refEl));
      for(int j = 0; j < refEl->getNumIPs(); j++){
        double dV = detJacs[j]*ipWeights->at(j);
        for(int k = 0; k < cell.size(); k++){
          for(int l = 0; l < cell.size(); l++){
            integral += ipShapes->at(j)[l]*(ipShapes->at(j)[k])*(field1->getValues()->at(i*(refEl->getNumNodes()) + k))*(field2->getValues()->at(i*(refEl->getNumNodes()) + l))*dV;
          }
        }
      }
    }
    return integral;
  };

  static std::vector< std::vector<double> > getRefNormals(int dim){
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

  static double analyticalAdvectionSolution(std::function<double(std::vector<double>)> initialSol, const std::vector<double> & velocity, const double t, const std::vector<double> & x){
    std::vector<double> effectiveX(x.size(), 0.0);
    EMap<EVector>(effectiveX.data(), effectiveX.size()) = EMap<const EVector>(x.data(), x.size()) - t*EMap<const EVector>(velocity.data(), velocity.size());
    return initialSol(effectiveX);
  };
  
  static double morelet(const std::vector<double> & x, const std::vector<double> & c, double freq, double dev){
    double r = (EMap<const EVector>(x.data(), x.size()) - EMap<const EVector>(c.data(), c.size())).norm();
    return (0.5*std::exp(-0.5*std::pow((r/dev), 2.0))*std::cos(freq*r));
  };

  static double moreletGrad(const std::vector<double> & x, const std::vector<double> & c, double freq, double dev, int k){
    double r = (EMap<const EVector>(x.data(), x.size()) - EMap<const EVector>(c.data(), c.size())).norm();
    return (0.5*((x[k]-c[k])/r)*std::exp(-0.5*std::pow((r/dev), 2.0))*((r/dev)*std::cos(freq*r) + freq*std::sin(freq*r)));
  };
  
  static std::vector<EVector> calculateOutwardNormal(Mesh * myMesh, int faceInd){
    const ReferenceElement * refEl = myMesh->getReferenceElement();
    const ReferenceElement * fEl = refEl->getFaceElement();
    int dimSpace = myMesh->getNodeSpaceDimension();
    int nNodes = fEl->getNumNodes();
    std::vector<int> faceIndexes;
    std::vector< std::vector<double> > faceNodes;
    const std::vector< std::vector< std::vector<double> > > * derivShapes = fEl->getDerivShapeFunctions();
    myMesh->getFace(faceInd, &faceIndexes);
    myMesh->getSlicePoints(faceIndexes, &faceNodes);
    EVector testVec = EMap<EVector>(faceNodes[0].data(), faceNodes[0].size());
    std::vector<double> node(dimSpace, 0.0);
    std::vector<int> face2Cell(2, 0);
    std::vector<int> cell(refEl->getNumNodes(), 0);
    myMesh->getFace2Cell(faceInd, &face2Cell);
    myMesh->getCell(face2Cell[0], &cell);
    std::vector<int>::iterator it;
    int innerNode = 0;
    for(int i = 0; i < cell.size(); i++){
      it = std::find(faceIndexes.begin(), faceIndexes.end(), cell[i]);
      if(it == faceIndexes.end()){
        innerNode = i;
        break;
      }
    }
    myMesh->getPoint(cell[innerNode], &node);
    testVec -= EMap<EVector>(node.data(), node.size());
    std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
    EMatrix jacobian(dimSpace, dimSpace - 1);
    for(int i = 0; i < nNodes; i++){
      jacobian *= 0.0;
      for(int j = 0; j < nNodes; j++){
        jacobian += EMap<const EVector>(faceNodes[j].data(), faceNodes[j].size())*(EMap<const EVector>(derivShapes->at(i)[j].data(), derivShapes->at(i)[j].size())).transpose();
      }
      res[i] = (jacobian*jacobian.transpose()).fullPivLu().kernel();
      res[i].normalize();
      if(res[i].dot(testVec) < 0){
        res[i] *= -1;
      }
    }
    return res;
  };


};//TestUtils


struct SimRun{
  std::string dim;
  std::string meshSize;
  std::string order;
  std::string timeStep;
  std::string meshLocation;
  std::string rk;
  double linAlgErr = 0;
  double l2Err = 0;
  double dL2Err = 0;
  std::chrono::duration<double> runtime = std::chrono::duration<double>::zero();
  std::chrono::duration<double> setup= std::chrono::duration<double>::zero();
  std::chrono::duration<double> assembly= std::chrono::duration<double>::zero();
  std::chrono::duration<double> resolution= std::chrono::duration<double>::zero();
  std::chrono::duration<double> post= std::chrono::duration<double>::zero();
};

}

#endif//TESTUTILS_H

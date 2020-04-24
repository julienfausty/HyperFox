#include <catch2/catch.hpp>
#include <string>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include "Mesh.h"
#include "Field.h"
#include "HDF5Io.h"
#include "Operator.h"
#include "LaplaceModel.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "CGSolver.h"
#include "TestUtils.h"

using namespace hfox;

struct SimRun{
  std::string dim;
  std::string meshSize;
  std::string order;
  std::string meshLocation;
  double linAlgErr;
  double l2Err;
  double dL2Err;
  std::chrono::duration<double> runtime;
};

double analyticalSolution(const std::vector<double> & point){
  if(point.size() == 2){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else if(point.size() == 3){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else{
    return 0;
  }
};

double l2ProjectionNodeField(Field * field1, Field * field2, Mesh * mesh){
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

void runSimulation(SimRun * thisRun){
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order + ".h5";
  thisRun->meshLocation += meshName; 
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
  int nNodesPerFace = myMesh.getReferenceElement()->getFaceElement()->getNumNodes();
  Field dirichlet(&myMesh, Face, nNodesPerFace, 1);
  //set dirichlet field to analytical solution
  const std::set<int> * boundary = myMesh.getBoundaryFaces();
  std::vector<int> cell(nNodesPerFace, 0);
  std::vector< std::vector<double> > nodes(nNodesPerFace, std::vector<double>(myMesh.getNodeSpaceDimension(), 0.0));
  for(auto it = boundary->begin(); it != boundary->end(); it++){
    myMesh.getFace(*it, &cell);
    myMesh.getSlicePoints(cell, &nodes);
    for(int j = 0; j < nNodesPerFace; j++){
      dirichlet.getValues()->at((*it)*nNodesPerFace+j) = analyticalSolution(nodes[j]);
    }
  }
  Field sol(&myMesh, Node, 1, 1);
  Field anaSol(&myMesh, Node, 1, 1);
  //calculate analytical solution
  std::vector<double> node(myMesh.getNodeSpaceDimension());
  for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
    myMesh.getPoint(iNode, &node);
    anaSol.getValues()->at(iNode) = analyticalSolution(node);
  }
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["Dirichlet"] = &dirichlet;
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  LaplaceModel lapMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.maxits = 20000;
  myOpts.rtol = 1e-16;
  myOpts.verbose = false;
  PetscInterface petsciface(myOpts);
  CGSolver mySolver;
  mySolver.setVerbosity(0);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&lapMod);
  mySolver.setBoundaryModel(&dirMod);
  mySolver.initialize();
  mySolver.allocate();
  mySolver.assemble();
  mySolver.solve();
  //get linalg err
  const KSP * ksp = petsciface.getKSP();
  KSPGetResidualNorm(*ksp, &(thisRun->linAlgErr));
  //calculate l2 err
  Field residual(&myMesh, Node, 1, 1);
  double sumRes = 0, sumAna = 0;
  for(int i = 0; i < myMesh.getNumberPoints(); i++){
    residual.getValues()->at(i) = anaSol.getValues()->at(i) - sol.getValues()->at(i);
    sumRes += std::pow(residual.getValues()->at(i), 2);
    sumAna += std::pow(anaSol.getValues()->at(i), 2);
  }
  thisRun->l2Err = std::sqrt(l2ProjectionNodeField(&residual, &residual, &myMesh));
  thisRun->dL2Err = std::sqrt(sumRes/sumAna);
};

TEST_CASE("Testing regression CGLaplace", "[regression][CG][Laplace]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  meshSizes["3"] = {"3e-1"};
  meshSizes["2"] = {"3e-1", "2e-1", "1e-1"};
  std::vector<std::string> orders = {"1", "2", "3"};
  std::vector<SimRun> simRuns;
  for(auto it = meshSizes.begin(); it != meshSizes.end(); it++){
    for(auto itMs = it->second.begin(); itMs != it->second.end(); itMs++){
      for(int o = 0; o < orders.size(); o++){
        SimRun thisRun;
        thisRun.dim = it->first;
        thisRun.meshSize = *itMs;
        thisRun.order = orders[o];
        simRuns.push_back(thisRun);
      }
    }
  }

  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runSimulation(&(*it));
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    CHECK(it->l2Err < 1e-2);
  }

};

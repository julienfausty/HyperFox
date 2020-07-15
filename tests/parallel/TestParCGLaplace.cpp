#include <catch2/catch.hpp>
#include <string>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include "Mesh.h"
#include "ZoltanPartitioner.h"
#include "Field.h"
#include "HDF5Io.h"
#include "Operator.h"
#include "LaplaceModel.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "CGSolver.h"
#include "TestUtils.h"

using namespace hfox;

double analyticalSolution(const std::vector<double> & point){
  if(point.size() == 2){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else if(point.size() == 3){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else{
    return 0;
  }
};

void runSimulation(SimRun * thisRun){
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order + ".h5";
  thisRun->meshLocation += meshName; 
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
  ZoltanPartitioner zPart(&myMesh);
  zPart.initialize();
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
  Field residual(&myMesh, Node, 1, 1);
  std::vector<Field*> fieldList = {&sol, &anaSol, &dirichlet, &residual};
  //calculate analytical solution
  std::vector<double> node(myMesh.getNodeSpaceDimension());
  for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
    myMesh.getPoint(iNode, &node);
    anaSol.getValues()->at(iNode) = analyticalSolution(node);
  }
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["Dirichlet"] = &dirichlet;
  zPart.setFields(fieldList);
  zPart.computePartition();
  zPart.update();
  std::cout << "out of partitioning" << std::endl;
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  LaplaceModel lapMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.maxits = 20000;
  myOpts.rtol = 1e-16;
  myOpts.verbose = false;
  PetscInterface petsciface(myOpts);
  CGSolver mySolver;
  mySolver.setVerbosity(false);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&lapMod);
  mySolver.setBoundaryModel(&dirMod);
  std::cout << "finished set up" << std::endl;
  mySolver.initialize();
  std::cout << "finished initialize" << std::endl;
  mySolver.allocate();
  std::cout << "finished allocation" << std::endl;
  mySolver.assemble();
  std::cout << "finished assembly" << std::endl;
  mySolver.solve();
  std::cout << "finished solve" << std::endl;
  zPart.updateSharedInformation();
  //get linalg err
  const KSP * ksp = petsciface.getKSP();
  KSPGetResidualNorm(*ksp, &(thisRun->linAlgErr));
  //calculate l2 err
  double sumRes = 0, sumAna = 0;
  for(int i = 0; i < myMesh.getNumberPoints(); i++){
    residual.getValues()->at(i) = anaSol.getValues()->at(i) - sol.getValues()->at(i);
    sumRes += std::pow(residual.getValues()->at(i), 2);
    sumAna += std::pow(anaSol.getValues()->at(i), 2);
  }
  thisRun->l2Err = std::sqrt(TestUtils::l2ProjectionNodeField(&residual, &residual, &myMesh));
  thisRun->dL2Err = std::sqrt(sumRes/sumAna);
  double errBuff = 0;
  MPI_Allreduce(&(thisRun->dL2Err), &errBuff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  thisRun->dL2Err = errBuff;
  //hdfio.setField("Solution", &sol);
  //std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/LaplaceConvergence/CG/";
  //hdfio.write(writePath + meshName);
};

TEST_CASE("Testing regression CGLaplace", "[parallel][CG][Laplace]"){
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
  //std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/LaplaceConvergence/";
  //std::string writeFile = "CG/CGLaplace.csv";
  //std::ofstream f; f.open(writePath + writeFile);
  //f << "dim, order, h, linAlgErr, l2Err, dL2Err, runtime\n";
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runSimulation(&(*it));
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    CHECK(it->l2Err < 1e-2);
    double timeBuff, runtime;
    timeBuff = it->runtime.count();
    MPI_Allreduce(&timeBuff, &runtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //f << it->dim << ", ";
    //f << it->order << ", ";
    //f << it->meshSize << ", ";
    //f << it->linAlgErr << ", ";
    //f << it->l2Err << ", ";
    //f << it->dL2Err << ", ";
    //f << it->runtime.count() << "\n";
  }
  //f.close();

};

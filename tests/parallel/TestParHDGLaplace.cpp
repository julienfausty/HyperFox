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
#include "HDGLaplaceModel.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "TestUtils.h"

using namespace hfox;

double sinExp(const std::vector<double> & point){
  if(point.size() == 2){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else if(point.size() == 3){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else{
    return 0;
  }
};

void runHDGSimulation(SimRun * thisRun){
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
      dirichlet.getValues()->at((*it)*nNodesPerFace+j) = sinExp(nodes[j]);
    }
  }
  int nNodes = myMesh.getReferenceElement()->getNumNodes();
  Field sol(&myMesh, Cell, nNodes, 1);
  Field flux(&myMesh, Cell, nNodes, myMesh.getNodeSpaceDimension());
  Field trace(&myMesh, Face, nNodesPerFace, 1);
  Field tau(&myMesh, Face, nNodesPerFace, 1);
  std::fill(tau.getValues()->begin(), tau.getValues()->end(), 1.0);
  Field anaSol(&myMesh, Cell, nNodes, 1);
  Field residual(&myMesh, Cell, nNodes, 1);
  std::vector<Field*> fieldList = {&sol, &flux, &trace, &tau, &anaSol, &dirichlet, &residual};
  //calculate analytical solution
  std::vector<double> node(myMesh.getNodeSpaceDimension());
  for(int iCell = 0; iCell < myMesh.getNumberCells(); iCell++){
    myMesh.getCell(iCell, &cell);
    for(int k = 0; k < cell.size(); k++){
      int iNode = zPart.global2LocalNode(cell[k]);
      if(iNode != -1){
        myMesh.getPoint(iNode, &node);
      } else {
        myMesh.getGhostPoint(cell[k], &node);
      }
      anaSol.getValues()->at(iCell*nNodes + k) = sinExp(node);
    }
  }
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["Flux"] = &flux;
  fieldMap["Trace"] = &trace;
  fieldMap["Tau"] = &tau;
  fieldMap["Dirichlet"] = &dirichlet;
  zPart.setFields(fieldList);
  zPart.computePartition();
  zPart.update();
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGLaplaceModel lapMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.maxits = 20000;
  myOpts.rtol = 1e-16;
  myOpts.verbose = false;
  PetscInterface petsciface(myOpts);
  HDGSolver mySolver;
  mySolver.setVerbosity(false);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&lapMod);
  mySolver.setBoundaryModel(&dirMod);
  mySolver.initialize();
  mySolver.allocate();
  mySolver.assemble();
  mySolver.solve();
  zPart.updateSharedInformation();
  //get linalg err
  const KSP * ksp = petsciface.getKSP();
  KSPGetResidualNorm(*ksp, &(thisRun->linAlgErr));
  //calculate l2 err
  cell.resize(nNodes);
  double sumRes = 0, sumAna = 0;
  for(int i = 0; i < myMesh.getNumberCells(); i++){
    for(int j = 0; j < nNodes; j++){
      residual.getValues()->at(i*nNodes + j) = anaSol.getValues()->at(i*nNodes + j) - sol.getValues()->at(i*nNodes + j);
      sumRes += std::pow(residual.getValues()->at(i*nNodes + j), 2);
      sumAna += std::pow(anaSol.getValues()->at(i*nNodes + j), 2);
    }
  }
  thisRun->l2Err = std::sqrt(TestUtils::l2ProjectionCellField(&residual, &residual, &myMesh));
  thisRun->dL2Err = std::sqrt(sumRes/sumAna);
  double errBuff = 0;
  MPI_Allreduce(&(thisRun->dL2Err), &errBuff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  thisRun->dL2Err = errBuff;
  hdfio.setField("Solution", &sol);
  hdfio.setField("Analytical", &anaSol);
  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/parallel/HDG/";
  //hdfio.write(writePath + meshName);
};

TEST_CASE("Testing regression HDGLaplace", "[parallel][HDG][Laplace]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"2e-1"};
  meshSizes["2"] = {"1e-1", "7e-2", "5e-2"};
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

  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/parallel/HDG/";
  std::string writeFile = "Laplace.csv";
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //std::ofstream f;
  //if(rank == 0){
    //f.open(writePath + writeFile);
    //f << "dim, order, h, linAlgErr, l2Err, dL2Err, runtime\n";
  //}
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGSimulation(&(*it));
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    CHECK(it->l2Err < 1e-2);
    double timeBuff, runtime;
    timeBuff = it->runtime.count();
    MPI_Allreduce(&timeBuff, &runtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //if(rank == 0){
      //f << it->dim << ", ";
      //f << it->order << ", ";
      //f << it->meshSize << ", ";
      //f << it->linAlgErr << ", ";
      //f << it->l2Err << ", ";
      //f << it->dL2Err << ", ";
      //f << runtime << "\n";
      //f << std::flush;
    //}
  }
  //if(rank == 0){
    //f << std::endl;
    //f.close();
  //}

};

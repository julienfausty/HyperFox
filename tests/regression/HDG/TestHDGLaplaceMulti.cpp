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
#include "HDGLaplaceModel.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "TestUtils.h"

using namespace hfox;

double sinExpMulti(const std::vector<double> & point){
  if(point.size() == 2){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else if(point.size() == 3){
    return (std::sin(point[0]) * std::exp(point[1]));
  } else{
    return 0;
  }
};

void runHDGLaplaceMulti(SimRun * thisRun){
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order + ".h5";
  thisRun->meshLocation += meshName; 
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
  int nNodesPerFace = myMesh.getReferenceElement()->getFaceElement()->getNumNodes();
  int dim =  myMesh.getNodeSpaceDimension();
  Field dirichlet(&myMesh, Face, nNodesPerFace, dim);
  //set dirichlet field to analytical solution
  const std::set<int> * boundary = myMesh.getBoundaryFaces();
  std::vector<int> cell(nNodesPerFace, 0);
  std::vector< std::vector<double> > nodes(nNodesPerFace, std::vector<double>(myMesh.getNodeSpaceDimension(), 0.0));
  for(auto it = boundary->begin(); it != boundary->end(); it++){
    myMesh.getFace(*it, &cell);
    myMesh.getSlicePoints(cell, &nodes);
    for(int j = 0; j < nNodesPerFace; j++){
      for(int k = 0; k < dim; k++){
        dirichlet.getValues()->at(((*it)*nNodesPerFace+j)*dim + k) = sinExpMulti(nodes[j]);
      }
    }
  }
  int nNodes = myMesh.getReferenceElement()->getNumNodes();
  Field sol(&myMesh, Cell, nNodes, dim);
  Field flux(&myMesh, Cell, nNodes, dim*myMesh.getNodeSpaceDimension());
  Field trace(&myMesh, Face, nNodesPerFace, dim);
  Field tau(&myMesh, Face, nNodesPerFace, dim*dim);
  for(int i = 0; i < myMesh.getNumberFaces(); i++){
    for(int j = 0; j < nNodesPerFace; j++){
      EMap<EMatrix>(tau.getValues()->data() + (i*nNodesPerFace + j)*dim*dim, dim, dim) = EMatrix::Identity(dim, dim);
    }
  }
  Field anaSol(&myMesh, Node, 1, dim);
  //calculate analytical solution
  std::vector<double> node(myMesh.getNodeSpaceDimension());
  for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
    myMesh.getPoint(iNode, &node);
    for(int k = 0; k < dim; k++){
      anaSol.getValues()->at(iNode*dim + k) = sinExpMulti(node);
    }
  }
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["Flux"] = &flux;
  fieldMap["Trace"] = &trace;
  fieldMap["Tau"] = &tau;
  fieldMap["Dirichlet"] = &dirichlet;
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
  //get linalg err
  const KSP * ksp = petsciface.getKSP();
  KSPGetResidualNorm(*ksp, &(thisRun->linAlgErr));
  //calculate l2 err
  Field residual(&myMesh, Cell, nNodes, 1);
  cell.resize(nNodes);
  double sumRes = 0, sumAna = 0;
  for(int i = 0; i < myMesh.getNumberCells(); i++){
    myMesh.getCell(i, &cell);
    for(int j = 0; j < nNodes; j++){
      for(int k = 0; k < dim; k++){
        residual.getValues()->at(i*nNodes + j) += std::pow(anaSol.getValues()->at(cell[j]*dim + k) - sol.getValues()->at((i*nNodes + j)*dim + k), 2);
        sumAna += std::abs(anaSol.getValues()->at(cell[j]*dim + k));
      }
      residual.getValues()->at(i*nNodes + j) = std::sqrt(residual.getValues()->at(i*nNodes + j));
      sumRes += residual.getValues()->at(i*nNodes + j);
    }
  }
  thisRun->l2Err = std::sqrt(TestUtils::l2ProjectionCellField(&residual, &residual, &myMesh));
  thisRun->dL2Err = std::sqrt(sumRes/sumAna);
  //hdfio.setField("Solution", &sol);
  //hdfio.setField("Analytical", &anaSol);
  //hdfio.setField("Residual", &residual);
  //std::string writePath = "/home/jfausty/workspace/Postprocess/results/LaplaceMulti/HDG/";
  //hdfio.write(writePath + meshName);
};

TEST_CASE("Testing regression HDGLaplace with multiple DOFs", "[regression][HDG][LaplaceMulti]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  meshSizes["3"] = {"2e-1"};
  meshSizes["2"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"1e-1"};
  std::vector<std::string> orders = {"1", "2", "3"};
  //std::vector<std::string> orders = {"3"};
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

  //std::string writePath = "/home/jfausty/workspace/Postprocess/results/LaplaceMulti/";
  //std::string writeFile = "HDG/HDGLaplace.csv";
  //std::ofstream f; f.open(writePath + writeFile);
  //f << "dim, order, h, linAlgErr, l2Err, dL2Err, runtime\n";
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGLaplaceMulti(&(*it));
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    CHECK(it->l2Err < 1e-1);
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

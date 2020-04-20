#include <catch2/catch.hpp>
#include <string>
#include <chrono>
#include <cmath>
#include "Mesh.h"
#include "Field.h"
#include "HDF5Io.h"
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
  std::chrono::duration<double> runtime;
};

double analyticalSolution(const std::vector<double> & point){
  if(point.size() == 2){
    return (std::cos(point[0]) * std::cosh(point[1]));
  } else if(point.size() == 3){
    return (std::cos(point[0]) * std::cosh(point[1]) * std::cos(point[2]));
  } else{
    return 0;
  }
};

void runSimulation(SimRun * thisRun){
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  thisRun->meshLocation += "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  thisRun->meshLocation += "_ord-" + thisRun->order + ".h5";
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
  Field dirichlet(&myMesh, Face, 1, 1);
  //set dirichlet field to analytical solution
  Field sol(&myMesh, Node, 1, 1);
  Field anaSol(&myMesh, Node, 1, 1);
  //calculate analytical solution
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["Dirichlet"] = &dirichlet;
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  LaplaceModel lapMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.verbose = true;
  PetscInterface petsciface(myOpts);
  CGSolver mySolver;
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
  //calculate l2 err
};

TEST_CASE("Testing regression CGLaplace", "[regression][CG][Laplace]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["2"] = {"1e-1", "7e-2", "5e-2", "2e-2", "1e-2"};
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1", "8e-2"};
  meshSizes["2"] = {"1e-1", "7e-2", "5e-2"};
  meshSizes["3"] = {"3e-1", "2e-1"};
  std::vector<std::string> orders = {"1", "2", "3", "4"};
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
    std::cout << it->runtime.count() << "s" << std::endl;
  }
};

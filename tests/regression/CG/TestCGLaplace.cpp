#include <catch2/catch.hpp>
#include <string>
#include <chrono>
#include "CGSolver.h"
#include "LaplaceModel.h"
#include "DirichletModel.h"
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

void runSimulation(SimRun * thisRun){
};

double analyticalSolution(const std::vector<double> & point){
};

TEST_CASE("Testing regression CGLaplace", "[regression][CG][Laplace]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  meshSizes["2"] = {"1e-1", "7e-2", "5e-2", "2e-2", "1e-2"};
  meshSizes["3"] = {"3e-1", "2e-1", "1e-1", "8e-2"};
  std::vector<std::string> orders = {"1", "2", "3", "4", "5"};
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
  }
};

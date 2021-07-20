#include <catch2/catch.hpp>
#include <string>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <boost/filesystem.hpp>
#include <fstream>
#include "Mesh.h"
#include "Field.h"
#include "ZoltanPartitioner.h"
#include "HDF5Io.h"
#include "HDGConnectionLaplacianModel.h"
#include "RungeKutta.h"
#include "IntegratedDirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "TestUtils.h"

using namespace hfox;

TEST_CASE("Testing regression cases for the HDGConnectionLaplacianModel", "[regression][HDG][ConnectionLaplacianModel]"){
  std::vector<std::string> hs = {"1e-0", "7e-1", "5e-1"};
  std::vector<std::string> orders = {"1", "2", "3", "4", "5"};
  for(int iH = 0; iH < hs.size(); iH++){
    for(int iO = 0; iO < orders.size(); iO){
    }
  }
};

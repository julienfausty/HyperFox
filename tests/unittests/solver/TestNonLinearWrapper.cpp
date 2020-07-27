#include <catch2/catch.hpp>
#include <string>
#include "NonLinearWrapper.h"
#include "CGSolver.h"
#include "Mesh.h"
#include "Field.h"
#include "HDF5Io.h"
#include "TestUtils.h"

using namespace hfox;

void linearizedQuadraticSolver(Field* currentF, Field* prevF, Solver * solver){
  for(int i = 0; i < currentF->getValues()->size(); i++){
    currentF->getValues()->at(i) = prevF->getValues()->at(i) - std::pow(prevF->getValues()->at(i), 2)/(2.0*(prevF->getValues()->at(i)));
  }
};

TEST_CASE("Testing the NonLinearWrapper", "[unit][solver][NonLinearWrapper]"){
  HDF5Io hdfio;
  Mesh m(2, 1, "simplex");
  hdfio.setMesh(&m);
  hdfio.load(TestUtils::getRessourcePath() + std::string("/meshes/lightTri.h5"));
  std::map<std::string, Field* > fieldMap;
  Field sol(&m, Node, 1, 1);
  Field interSol(&m, Node, 1, 1);
  CGSolver cgSolve;
  std::vector<double> startingVals = {1.0, 2.0, -1.0, 0.25, 42.0, 1e-8};
  CHECK_NOTHROW(NonLinearWrapper());
  NonLinearWrapper wrap;
  CHECK_NOTHROW(wrap.setVerbosity(0));
  CHECK_NOTHROW(wrap.setSolutionFields(&sol, &interSol));
  CHECK_NOTHROW(wrap.setSolver(&cgSolve));
  CHECK_NOTHROW(wrap.setLinearizedSolver([&sol, &interSol](Solver * solver){linearizedQuadraticSolver(&sol, &interSol, solver);}));
  for(int i = 0; i < startingVals.size(); i++){
    SECTION("Testing starting value " + std::to_string(startingVals[i])){
      for(int k = 0; k < interSol.getValues()->size(); k++){
        interSol.getValues()->at(k) = startingVals[i];
        sol.getValues()->at(k) = 0.0;
      }
      CHECK_NOTHROW(wrap.solve());
      CHECK(wrap.getResidual() < 1e-6);
      for(int k = 0; k < interSol.getValues()->size(); k++){
        CHECK(std::abs(interSol.getValues()->at(k)) < 1e-4);
      }
    };
  }
};

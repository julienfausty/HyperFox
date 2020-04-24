#include <catch2/catch.hpp>
#include <string>
#include "CGSolver.h"
#include "Mesh.h"
#include "Field.h"
#include "DirichletModel.h"
#include "LaplaceModel.h"
#include "PetscInterface.h"
#include "HDF5Io.h"
#include "TestUtils.h"

using namespace hfox;

TEST_CASE("Testing the CGSolver", "[unit][solver][CGSolver]"){

  HDF5Io hdfio;
  Mesh m(2, 1, "simplex");
  hdfio.setMesh(&m);
  hdfio.load(TestUtils::getRessourcePath() + std::string("/meshes/lightTri.h5"));
  std::map<std::string, Field* > fieldMap;
  Field sol(&m, Node, 1, 1);
  fieldMap["Solution"] = &sol;
  Field dir(&m, Face, m.getReferenceElement()->getFaceElement()->getNumNodes(), 1);
  fieldMap["Dirichlet"] = &dir;
  std::fill(dir.getValues()->begin(), dir.getValues()->end(), 3.0);
  DirichletModel dirMod(m.getReferenceElement()->getFaceElement());
  LaplaceModel lapMod(m.getReferenceElement());
  PetscOpts myOpts;
  myOpts.verbose = false;
  PetscInterface petscIFace(myOpts);
  SECTION("Testing setup"){
    CGSolver cgSolve;

    CHECK_NOTHROW(cgSolve.setVerbosity(0));

    CHECK_THROWS(cgSolve.solve());
    CHECK_THROWS(cgSolve.assemble());
    CHECK_THROWS(cgSolve.allocate());

    CHECK_NOTHROW(cgSolve.setMesh(&m));
    CHECK_THROWS(cgSolve.solve());
    CHECK_THROWS(cgSolve.assemble());
    CHECK_THROWS(cgSolve.allocate());

    CHECK_NOTHROW(cgSolve.setFieldMap(&fieldMap));
    CHECK_THROWS(cgSolve.solve());
    CHECK_THROWS(cgSolve.assemble());
    CHECK_THROWS(cgSolve.allocate());

    CHECK_NOTHROW(cgSolve.setLinSystem(&petscIFace));
    CHECK_THROWS(cgSolve.solve());
    CHECK_THROWS(cgSolve.assemble());
    CHECK_THROWS(cgSolve.allocate());

    CHECK_NOTHROW(cgSolve.setModel(&lapMod));
    CHECK_THROWS(cgSolve.solve());
    CHECK_THROWS(cgSolve.assemble());
    CHECK_THROWS(cgSolve.allocate());

    CHECK_NOTHROW(cgSolve.setBoundaryModel(&dirMod));
    CHECK_THROWS(cgSolve.solve());
    CHECK_THROWS(cgSolve.assemble());

    CHECK_THROWS(cgSolve.allocate());
    CHECK_THROWS(cgSolve.solve());

    CHECK_NOTHROW(cgSolve.initialize());
    CHECK_NOTHROW(cgSolve.allocate());
    CHECK_THROWS(cgSolve.solve());

    CHECK_NOTHROW(cgSolve.assemble());
    CHECK_NOTHROW(cgSolve.solve());

    const std::vector<double> * vals = sol.getValues();
    for(int i = 0; i < m.getNumberPoints(); i++){
      CHECK((*vals)[i] == Approx(3.0).margin(1e-12));
    }
  };
};

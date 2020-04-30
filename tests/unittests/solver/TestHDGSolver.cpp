#include <catch2/catch.hpp>
#include <string>
#include "HDGSolver.h"
#include "Mesh.h"
#include "Field.h"
#include "DirichletModel.h"
#include "HDGLaplaceModel.h"
#include "PetscInterface.h"
#include "HDF5Io.h"
#include "TestUtils.h"

using namespace hfox;

TEST_CASE("Testing the HDGSolver", "[unit][solver][HDGSolver]"){

  HDF5Io hdfio;
  Mesh m(2, 1, "simplex");
  hdfio.setMesh(&m);
  hdfio.load(TestUtils::getRessourcePath() + std::string("/meshes/lightTri.h5"));
  std::map<std::string, Field* > fieldMap;
  Field sol(&m, Cell, m.getReferenceElement()->getNumNodes(), 1);
  Field flux(&m, Cell, m.getReferenceElement()->getNumNodes(), 2);
  Field dir(&m, Face, m.getReferenceElement()->getFaceElement()->getNumNodes(), 1);
  Field lambda(&m, Face, m.getReferenceElement()->getFaceElement()->getNumNodes(), 1);
  std::fill(dir.getValues()->begin(), dir.getValues()->end(), 3.0);
  DirichletModel dirMod(m.getReferenceElement()->getFaceElement());
  HDGLaplaceModel hdgLapMod(m.getReferenceElement());
  PetscOpts myOpts;
  myOpts.verbose = true;
  PetscInterface petscIFace(myOpts);
  SECTION("Testing setup"){
    HDGSolver hdgSolve;

    CHECK_NOTHROW(hdgSolve.setVerbosity(1));

    CHECK_THROWS(hdgSolve.solve());
    CHECK_THROWS(hdgSolve.assemble());
    CHECK_THROWS(hdgSolve.allocate());

    CHECK_NOTHROW(hdgSolve.setMesh(&m));
    CHECK_THROWS(hdgSolve.solve());
    CHECK_THROWS(hdgSolve.assemble());
    CHECK_THROWS(hdgSolve.allocate());

    fieldMap["Solution"] = &sol;
    fieldMap["Flux"] = &flux;
    fieldMap["Trace"] = &lambda;
    fieldMap["Dirichlet"] = &dir;
    CHECK_NOTHROW(hdgSolve.setFieldMap(&fieldMap));
    CHECK_THROWS(hdgSolve.solve());
    CHECK_THROWS(hdgSolve.assemble());
    CHECK_THROWS(hdgSolve.allocate());

    CHECK_NOTHROW(hdgSolve.setLinSystem(&petscIFace));
    CHECK_THROWS(hdgSolve.solve());
    CHECK_THROWS(hdgSolve.assemble());
    CHECK_THROWS(hdgSolve.allocate());

    CHECK_NOTHROW(hdgSolve.setModel(&hdgLapMod));
    CHECK_THROWS(hdgSolve.solve());
    CHECK_THROWS(hdgSolve.assemble());
    CHECK_THROWS(hdgSolve.allocate());

    CHECK_NOTHROW(hdgSolve.setBoundaryModel(&dirMod));
    CHECK_THROWS(hdgSolve.solve());
    CHECK_THROWS(hdgSolve.assemble());

    CHECK_THROWS(hdgSolve.allocate());
    CHECK_THROWS(hdgSolve.solve());

    CHECK_NOTHROW(hdgSolve.initialize());
    CHECK_NOTHROW(hdgSolve.allocate());
    CHECK_THROWS(hdgSolve.solve());

    CHECK_NOTHROW(hdgSolve.assemble());
    CHECK_NOTHROW(hdgSolve.solve());

    const std::vector<double> * solVals, *fluxVals, *traceVals;
    solVals = sol.getValues();
    fluxVals = flux.getValues();
    traceVals = lambda.getValues();
    int nNodesPerEl = m.getReferenceElement()->getNumNodes();
    int nNodesPerFc = m.getReferenceElement()->getFaceElement()->getNumNodes();
    for(int i = 0; i < m.getNumberCells(); i++){
      for(int j = 0; j < nNodesPerEl; j++){
        CHECK((*solVals)[i*nNodesPerEl + j] == Approx(3.0).margin(1e-12));
        for(int k = 0; k < 2; k++){
          CHECK((*fluxVals)[i*nNodesPerEl + j*2 + k] == Approx(0.0).margin(1e-12));
        }
      }
    }
    for(int i = 0; i < m.getNumberFaces(); i++){
      for(int j = 0; j < nNodesPerFc; j++){
        CHECK((*traceVals)[i*nNodesPerFc + j] == Approx(3.0).margin(1e-12));
      }
    }
  };
};

#include <catch2/catch.hpp>
#include <string>
#include <numeric>
#include "DirichletModel.h"
#include "AssemblyType.h"


using namespace hfox;

TEST_CASE("Testing the DirichletModel", "[unit][model][DirichletModel]"){

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 5; j++){
      SECTION("Testing setup (dim = " + std::to_string(i+1) + ", order = " + std::to_string(j + 1) + ")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        DirichletModel dirModel(&refEl);
        CHECK_THROWS(dirModel.compute());
        std::map<std::string, std::vector<double> > fm;
        CHECK_THROWS(dirModel.setFieldMap(&fm));
        std::vector<double> dirichlet(refEl.getNumNodes(), 0.0);
        std::iota(dirichlet.begin(), dirichlet.end(), 0.0);
        fm["Dirichlet"] = dirichlet;
        CHECK_NOTHROW(dirModel.setFieldMap(&fm));
        CHECK_THROWS(dirModel.compute());
        CHECK_NOTHROW(dirModel.allocate(1));
        CHECK_NOTHROW(dirModel.compute());
        CHECK(dirModel.getAssemblyType()->matrix == Set);
        CHECK(dirModel.getAssemblyType()->rhs == Set);
      };
    }
  }

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 5; j++){
      SECTION("Testing local objects (dim = " + std::to_string(i+1) + ", order = " + std::to_string(j + 1) + ")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        DirichletModel dirModel(&refEl);
        dirModel.allocate(1);
        std::map<std::string, std::vector<double> > fm;
        std::vector<double> dirichlet(refEl.getNumNodes(), 0.0);
        std::iota(dirichlet.begin(), dirichlet.end(), 0.0);
        fm["Dirichlet"] = dirichlet;
        dirModel.setFieldMap(&fm);
        dirModel.compute();
        CHECK(*(dirModel.getLocalMatrix()) == EMatrix::Identity(refEl.getNumNodes(), refEl.getNumNodes()));
        CHECK(*(dirModel.getLocalRHS()) == Eigen::Map<EVector>(dirichlet.data(), dirichlet.size()));
      };
    }
  }
};


#include "catch.hpp"
#include "Simplex.h"

using namespace hfox;

TEST_CASE("Tests for Simplex class.", "[Simplex][CPolytope][Polytope][unit][mesh]"){
  std::vector< std::vector<double> > testpoints2;
  int i, j;
  for(i = 0; i < 3; i++){
    for(j = 0; j <2; j++){
      testpoints2[i][j] = i+j;
    }
  }
  std::vector< std::vector<double> > testpoints3;
  for(i = 0; i < 4; i++){
    for(j = 0; j <3; j++){
      testpoints2[i][j] = i+j;
    }
  }
  SECTION("Construction."){
    REQUIRE_NOTHROW(new Simplex<1>());
    REQUIRE_NOTHROW(new Simplex<2>());
    REQUIRE_NOTHROW(new Simplex<3>());
    REQUIRE_NOTHROW(new Simplex<4>());
  };
};

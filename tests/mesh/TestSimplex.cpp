#include "catch.hpp"
#include <iostream>
#include "Simplex.h"

using namespace hfox;

TEST_CASE("Tests for Simplex class.", "[Simplex][CPolytope][Polytope][unit][mesh]"){
  std::vector< std::vector<double> > testpoints2;
  testpoints2.resize(3);
  int i, j;
  for(i = 0; i < 3; i++){
    (testpoints2[i]).resize(2);
    for(j = 0; j <2; j++){
      (testpoints2[i])[j] = i+j;
    }
  }
  std::vector< std::vector<int> > testfaces2 = {{2,1},{0,2},{1,0}};
  std::vector< std::vector<double> > testpoints3;
  testpoints3.resize(4);
  for(i = 0; i < 4; i++){
    (testpoints3[i]).resize(3);
    for(j = 0; j <3; j++){
      (testpoints3[i])[j] = i+j;
    }
  }
  std::vector< std::vector<int> > testfaces3={{1,2,3},
    {3,2,0},{0,1,3},{2,1,0}};
  std::vector< std::vector<double> > testpoints4;
  testpoints4.resize(5);
  for(i = 0; i < 5; i++){
    (testpoints4[i]).resize(4);
    for(j = 0; j <4; j++){
      (testpoints4[i])[j] = i+j;
    }
  }
  std::vector< std::vector<int> > testfaces4={{4,3,2,1},
    {0,2,3,4},{4,3,1,0},{0,1,2,4},{3,2,1,0}};
  std::vector<std::vector<double> > badtestpoints2 = testpoints2;
  badtestpoints2.resize(2);
  std::vector<std::vector<double> > badtestpoints3 = testpoints3;
  for(i = 0; i < 4; i++){
    (badtestpoints3[i]).resize(2);
  }
  SECTION("Construction."){
    REQUIRE_NOTHROW(new Simplex<1>());
    REQUIRE_NOTHROW(new Simplex<2>());
    REQUIRE_NOTHROW(new Simplex<3>());
    REQUIRE_NOTHROW(new Simplex<4>());
    REQUIRE_NOTHROW(new Simplex<2>(testpoints2));
    REQUIRE_NOTHROW(new Simplex<3>(testpoints3));
    REQUIRE_NOTHROW(new Simplex<4>(testpoints4));
    Simplex<2> sim2(testpoints2);
    std::vector<std::vector<double> > r2 = *(sim2.getPoints());
    std::vector<std::vector<int> > i2 = *(sim2.getFaces());
    CHECK(sim2.getDimension() == 2);
    CHECK(r2 == testpoints2);
    CHECK(i2 == testfaces2);
    Simplex<3> sim3(testpoints3);
    std::vector<std::vector<double> > r3 = *(sim3.getPoints());
    std::vector<std::vector<int> > i3 = *(sim3.getFaces());
    CHECK(sim3.getDimension() == 3);
    CHECK(r3 == testpoints3);
    CHECK(i3 == testfaces3);
    Simplex<4> sim4(testpoints4);
    std::vector<std::vector<double> > r4 = *(sim4.getPoints());
    std::vector<std::vector<int> > i4 = *(sim4.getFaces());
    CHECK(sim4.getDimension() == 4);
    CHECK(r4 == testpoints4);
    CHECK(i4 == testfaces4);
    REQUIRE_THROWS(new Simplex<2>(badtestpoints2));
    REQUIRE_THROWS(new Simplex<3>(badtestpoints3));
    REQUIRE_THROWS(new Simplex<2>(badtestpoints3));
  };
};

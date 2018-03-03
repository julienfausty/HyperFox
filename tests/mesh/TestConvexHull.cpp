#include "catch.hpp"
#include "ConvexHull.h"

using namespace hfox;

TEST_CASE("Testing the ConvexHull class.", "[ConvexHull][tool][unit][mesh]"){
  //Database
  std::vector< std::vector<double> > testpoints1 = {{0,0}, {1,0}, {0,1}};
  std::vector< std::vector<int> > testfaces1 = {{2,1},{0,2},{1,0}};
  std::vector< std::vector<double> > testpoints2 = {{0,0}, {1,0}, {0,1},
    {0.25,0.25}};
  std::vector< std::vector<int> > testfaces2 = {{2,1},{0,2},{1,0}};
  std::vector< std::vector<double> > testpoints3 = {
    {1,0}, 
    {0.86602540378,0.5}, 
    {0.5,0.86602540378}, 
    {0, 1}, 
    {-0.5,0.86602540378},
    {-0.86602540378, 0.5},
    {-1, 0},
    {-0.86602540378, -0.5},
    {-0.5, -0.86602540378},
    {0, -1},
    {0.5, -0.86602540378},
    {0.86602540378, -0.5}
  };
  std::vector< std::vector<int> > testfaces3 = {
    {0,1},{1,2},{2,3},
    {3,4},{4,5},{5,6},
    {6,7},{7,8},{7,9},
    {9,10},{10,11},{11,12}
  };
  std::vector< std::vector<double> > testpoints4 = {{0,0,0}, {1,0,0}, {0,1,0},
    {0,0,1}};
  std::vector< std::vector<int> > testfaces4 = {{1,2,3},{3,2,0},{0,1,3},
    {2,1,0}};
  std::vector< std::vector<double> > testpoints5 = {{0,0,0}, {1,0,0}, {1,1,0}, 
    {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};
  // just test the number of faces for case 5 -> should be 8*2=16.
  std::vector< std::vector<double> > testpoints6 = {{0,0,0}, {1,0,0}, {0,1,0},
    {0,0,1},{1,1,1}};
  std::vector< std::vector<int> > testfaces6 = {
    {3,2,0},{0,1,3},{2,1,0},
    {2,3,4},{4,3,1},{1,2,4}};
    std::vector< std::vector<double> > * ppoints;
    std::vector< std::vector<int> > * pfaces;
  SECTION("Testing ConvexHull class construction."){
    REQUIRE_NOTHROW(new ConvexHull());
    ConvexHull ch;
    CHECK(ch.getAlgoType() == QuickHull);
  };
  SECTION("Testing Normal computation."){
    // pre test
    ConvexHull ch;
    std::vector<int> chosenface;
    EVector realnormal(2);
    EVector n(2);
    EVector diff(2);
    //tests
    ppoints = &testpoints1;
    ch.setVertexes(ppoints);
    chosenface.resize(2); chosenface[0] = 2; chosenface[1] = 1;
    realnormal.resize(2); realnormal(0) = 0.707107; realnormal(1) = 0.707107;
    n = ch.calcNormalVector(chosenface);
    diff = realnormal - n;
    CHECK(diff.dot(diff)<1e-12);
    ppoints = &testpoints2;
    ch.setVertexes(ppoints);
    chosenface.resize(2); chosenface[0] = 0; chosenface[1] = 2;
    realnormal.resize(2); realnormal(0) = -1; realnormal(1) = 0;
    n = ch.calcNormalVector(chosenface);
    diff = realnormal - n;
    CHECK(diff.dot(diff)<1e-12);
    ppoints = &testpoints3;
    ch.setVertexes(ppoints);
    chosenface.resize(2); chosenface[0] = 2; chosenface[1] = 1;
    realnormal.resize(2); realnormal(0) = 0.70710678118; 
    realnormal(1) = 0.70710678118;
    n = ch.calcNormalVector(chosenface);
    diff = realnormal - n;
    CHECK(diff.dot(diff)<1e-12);
    ppoints = &testpoints4;
    ch.setVertexes(ppoints);
    chosenface.resize(3); chosenface[0] = 1; chosenface[1] = 2; 
    chosenface[2] = 3;
    realnormal.resize(3); realnormal(0) = 0.57735026919; 
    realnormal(1) = 0.57735026919; realnormal(2) = 0.57735026919;
    n.resize(3); diff.resize(3);
    n = ch.calcNormalVector(chosenface);
    diff = realnormal - n;
    CHECK(diff.dot(diff)<1e-12);
    ppoints = &testpoints5;
    ch.setVertexes(ppoints);
    chosenface.resize(3); chosenface[0] = 1; chosenface[1] = 2; 
    chosenface[2] = 3;
    realnormal.resize(3); realnormal(0) = 0.0; 
    realnormal(1) = 0.0; realnormal(2) = 1;
    n.resize(3); diff.resize(3);
    n = ch.calcNormalVector(chosenface);
    diff = realnormal - n;
    CHECK(diff.dot(diff)<1e-12);
    ppoints = &testpoints6;
    ch.setVertexes(ppoints);
    chosenface.resize(3); chosenface[0] = 1; chosenface[1] = 2; 
    chosenface[2] = 4;
    realnormal.resize(3); realnormal(0) = 0.57735026919; 
    realnormal(1) = 0.57735026919; realnormal(2) = -0.57735026919;
    n.resize(3); diff.resize(3);
    n = ch.calcNormalVector(chosenface);
    diff = realnormal - n;
    CHECK(diff.dot(diff)<1e-12);
  };
  SECTION("Testing distance to face."){
    // pre test
    ConvexHull ch;
    std::vector<int> chosenface;
    EVector point(2);
    double distance;
    double diff;
  };
  SECTION("Testing Quickhull."){
    ConvexHull ch;
    ch.setFaces(pfaces);
    ppoints = &testpoints1;
    ch.setVertexes(ppoints);
    ch.computeConvexHull();
    CHECK(ch.getDimension() == 2);
    CHECK((*pfaces) == testfaces1);
    ppoints = &testpoints2;
    ch.setVertexes(ppoints);
    ch.computeConvexHull();
    CHECK(ch.getDimension() == 2);
    CHECK((*pfaces) == testfaces2);
    ppoints = &testpoints3;
    ch.setVertexes(ppoints);
    ch.computeConvexHull();
    CHECK(ch.getDimension() == 2);
    CHECK((*pfaces) == testfaces3);
    ppoints = &testpoints4;
    ch.setVertexes(ppoints);
    ch.computeConvexHull();
    CHECK(ch.getDimension() == 3);
    CHECK((*pfaces) == testfaces4);
    ppoints = &testpoints5;
    ch.setVertexes(ppoints);
    ch.computeConvexHull();
    CHECK(ch.getDimension() == 3);
    CHECK((*pfaces).size() == 16);
    ppoints = &testpoints6;
    ch.setVertexes(ppoints);
    ch.computeConvexHull();
    CHECK(ch.getDimension() == 3);
    CHECK((*pfaces) == testfaces6);
  };
};

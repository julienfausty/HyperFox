#include "catch.hpp"
#include <string>
#include <exception>
#include <vector>
#include "ErrorHandle.h"

using namespace hfox;

TEST_CASE("Tests for ErrorHandle class.", "[ErrorHandle][io][checkers][unit]"){
  std::string err_msg("Test error message.");
  std::string err_func("testErrFunc");
  std::string err_class("TestErrClass");

  SECTION("Test constructors"){
    CHECK_NOTHROW(new ErrorHandle());
    CHECK_NOTHROW(new ErrorHandle(err_msg));
    CHECK_NOTHROW(new ErrorHandle(err_class, err_func));
    CHECK_NOTHROW(new ErrorHandle(err_class, err_func, err_msg));
    CHECK_NOTHROW(new ErrorHandle("test msg"));
  };
  SECTION("Test what"){
    ErrorHandle errhand1(err_msg);
    ErrorHandle errhand2(err_class, err_func);
    ErrorHandle errhand3(err_class, err_func, err_msg);
    std::string delim = " : ";
    std::string temp;
    std::string place(errhand1.what());
    CHECK(place == ("SomeClass : someFunction : " + err_msg));
    temp = err_class + delim + err_func + delim + "something went wrong.\n";
    place = errhand2.what();
    CHECK(place == temp);
    temp = err_class + delim + err_func + delim + err_msg;
    place = errhand3.what();
    CHECK(place == temp);
  };
  SECTION("Test checkPointsList"){
    std::vector< std::vector<double> > goodtestpoints(2);
    std::vector<double> p1 = {1, 0.1};
    std::vector<double> p2 = {0.5, 10};
    std::vector< std::vector<double> > badtestpoints(2);
    std::vector<double> p2prime = {1, 3, 4};
    goodtestpoints[0] = p1;
    goodtestpoints[1] = p2;
    badtestpoints[0] = p1;
    badtestpoints[1] = p2prime;
    ErrorHandle errhand("TestClass", "testFunc");
    CHECK_NOTHROW(errhand.checkPointList(goodtestpoints));
    try{
      errhand.checkPointList(badtestpoints);
    }catch(ErrorHandle & eh){
      std::string w = eh.what();
      std::string r = "TestClass : testFunc : point list dimensions not respected! (2 != 3)\n";
      CHECK(w == r);
    }
    CHECK_THROWS(errhand.checkPointList(badtestpoints));
  };

  SECTION("Test checkIndexList"){
    std::vector< std::vector<int> > goodtestindexes(2);
    std::vector<int> p1 = {1, 2};
    std::vector<int> p2 = {3, 4};
    std::vector< std::vector<int> > badtestindexes(2);
    std::vector<int> p2prime = {1, 10};
    goodtestindexes[0] = p1;
    goodtestindexes[1] = p2;
    badtestindexes[0] = p1;
    badtestindexes[1] = p2prime;
    ErrorHandle errhand("TestClass", "testFunc");
    CHECK_NOTHROW(errhand.checkIndexList(9, goodtestindexes));
    try{
      errhand.checkIndexList(9, badtestindexes);
    }catch(ErrorHandle & eh){
      std::string w = eh.what();
      std::string r = "TestClass : testFunc : index list range not respected! (10 > 9)\n";
      CHECK(w == r);
    }
    CHECK_THROWS(errhand.checkIndexList(9, badtestindexes));
  };
};

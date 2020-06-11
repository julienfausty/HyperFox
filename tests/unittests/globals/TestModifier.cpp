#include <catch2/catch.hpp>
#include <vector>
#include <set>
#include <iostream>
#include "Modifier.h"

using namespace hfox;

TEST_CASE("Testing the Modifier class", "[unit][globals][Modifier]"){

  SECTION("Construction"){
    CHECK_NOTHROW(Modifier< std::vector<int> >());
    CHECK_NOTHROW(Modifier< std::vector<double> >());
  }

  int maxIt = 5;

  SECTION("Append"){
    Modifier< std::vector<int> > mod;
    mod.setType(APPEND);
    for(int i = 0; i < maxIt; i++){
      std::vector<int> vals(i);
      std::iota(vals.begin(), vals.end(), 0);
      mod.setValues(&vals);
      for(int j = 0; j < maxIt; j++){
        std::vector<int> v(j);
        std::iota(v.begin(), v.end(), 0);
        mod.modify(&v);
        CHECK(v.size() == i + j);
        for(int l = 0; l < j; l++){
          CHECK(v[l] == l);
        }
        for(int l = 0; l < i; l++){
          CHECK(v[j + l] == l);
        }
      }
    }
    Modifier< std::vector<double> > mod2;
    mod2.setType(APPEND); 
    for(int i = 0; i < maxIt; i++){
      std::vector<double> vals(i);
      std::iota(vals.begin(), vals.end(), 0.0);
      mod2.setValues(&vals);
      for(int j = 0; j < maxIt; j++){
        std::vector<double> v(j);
        std::iota(v.begin(), v.end(), 0);
        mod2.modify(&v);
        CHECK(v.size() == i + j);
        for(int l = 0; l < j; l++){
          CHECK(v[l] == l);
        }
        for(int l = 0; l < i; l++){
          CHECK(v[j + l] == l);
        }
      }
    }
  }

  SECTION("Remove"){
    Modifier< std::vector<int> > mod;
    mod.setType(REMOVE);
    for(int i = 0; i < maxIt; i++){
      std::vector<int> iv(100);
      std::iota(iv.begin(), iv.end(), 0);
      std::vector<double> dv(100);
      std::iota(dv.begin(), dv.end(), 0);
      std::vector<int> vals(10);
      std::iota(vals.begin(), vals.end(), i);
      mod.setValues(&vals);
      mod.modify(&iv);
      mod.modify(&dv);
      for(int k = 0; k < vals.size(); k++){
        CHECK(std::find(iv.begin(), iv.end(), vals[k]) == iv.end());
        CHECK(std::find(dv.begin(), dv.end(), vals[k]) == dv.end());
      }
    }
    for(int i = 0; i < maxIt; i++){
      std::vector<int> iv(100);
      std::iota(iv.begin(), iv.end(), 0);
      std::vector<double> dv(100);
      std::iota(dv.begin(), dv.end(), 0);
      std::vector<int> vals(10);
      std::iota(vals.begin(), vals.end(), i);
      for(int k = 0; k < vals.size(); k++){
        vals[k] *= 2;
      }
      mod.setValues(&vals);
      mod.modify(&iv);
      mod.modify(&dv);
      for(int k = 0; k < vals.size(); k++){
        CHECK(std::find(iv.begin(), iv.end(), vals[k]) == iv.end());
        CHECK(std::find(dv.begin(), dv.end(), vals[k]) == dv.end());
      }
    }
  };
};


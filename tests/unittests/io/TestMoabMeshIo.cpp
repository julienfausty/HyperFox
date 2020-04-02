#include <catch2/catch.hpp>
#include <iostream>
#include "configRessources.h"
#include "MoabMeshIo.h"

using namespace hfox;

TEST_CASE("Unit testing the MoabMeshIo class", "[unit][io][MoabMeshIo]"){
  std::string meshDirPath = getRessourcePath() + "/meshes/";
  std::cout << meshDirPath << std::endl;
};

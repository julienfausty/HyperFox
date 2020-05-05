#include <catch2/catch.hpp>
#include "HDGBase.h"
#include "HDGDiffusion.h"

using namespace hfox;

TEST_CASE("Testing HDGDiffusion operator", "[unit][operator][HDGDiffusion]"){

  SECTION("Testing allocation necessary"){
    ReferenceElement refEl(2, 1, "simplex");
    HDGBase base(&refEl);
    HDGDiffusion diff(&refEl);
    std::vector<double> dummydV(refEl.getNumIPs(), 0.0);
    std::vector<EMatrix> dummyJac(refEl.getNumIPs(), EMatrix::Identity(1,1));
    CHECK_THROWS(diff.assemble(dummydV, dummyJac));
    CHECK_NOTHROW(diff.allocate(1));
    CHECK_THROWS(diff.assemble(dummydV, dummyJac));
    base.allocate(1);
    base.calcNormals(*(refEl.getNodes()), dummyJac);
    CHECK_NOTHROW(diff.setFromBase(base.getNormals()));
    CHECK_NOTHROW(diff.assemble(dummydV, dummyJac));
  };
};

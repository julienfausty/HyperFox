#include <catch2/catch.hpp>
#include <string>
#include "LaplaceModel.h"
#include "AssemblyType.h"
#include "Operator.h"
#include "Diffusion.h"


using namespace hfox;

TEST_CASE("Testing the LaplaceModel", "[unit][model][LaplaceModel]"){

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 5; j++){
      SECTION("Testing setup(dim = " + std::to_string(i+1) + ", order = " + std::to_string(j + 1) + ")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        LaplaceModel lapModel(&refEl);
        CHECK_THROWS(lapModel.compute());
        std::map<std::string, std::vector<double> > fm;
        CHECK_NOTHROW(lapModel.setFieldMap(&fm));
        fm["DiffusionTensor"] = std::vector<double>(refEl.getNumNodes(), 0.0);
        CHECK_NOTHROW(lapModel.setFieldMap(&fm));
        CHECK_THROWS(lapModel.compute());
        CHECK_NOTHROW(lapModel.allocate(1));
        CHECK_THROWS(lapModel.compute());
        CHECK_NOTHROW(lapModel.setElementNodes(refEl.getNodes()));
        CHECK_NOTHROW(lapModel.compute());
        CHECK(lapModel.getAssemblyType()->matrix == Add);
        CHECK(lapModel.getAssemblyType()->rhs == None);
      };

      SECTION("Testing no diff coeff Matrix(dim = " + std::to_string(i+1) + ", order = " + std::to_string(j + 1) + ")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        Diffusion diffOp(&refEl);
        diffOp.allocate(1);
        std::vector<EMatrix> jacs = Operator::calcJacobians(*(refEl.getNodes()), &refEl);
        std::vector<EMatrix> invJacs = Operator::calcInvJacobians(jacs);
        std::vector<double> dV = Operator::calcMeasure(Operator::calcDetJacobians(jacs), &refEl);
        diffOp.assemble(dV, invJacs);
        LaplaceModel lapModel(&refEl);
        lapModel.allocate(1);
        lapModel.setElementNodes(refEl.getNodes());
        std::map<std::string, std::vector<double> > fm;
        lapModel.setFieldMap(&fm);
        lapModel.compute();
        CHECK(*(lapModel.getLocalMatrix()) == *(diffOp.getMatrix()));
        CHECK(*(lapModel.getLocalRHS()) == EVector::Zero(refEl.getNumNodes()));
      };

      SECTION("Testing with diff coeff Matrix(dim = " + std::to_string(i+1) + ", order = " + std::to_string(j + 1) + ")"){
        ReferenceElement refEl(i+1, j+1, "simplex");
        Diffusion diffOp(&refEl);
        diffOp.allocate(1);
        std::vector<EMatrix> jacs = Operator::calcJacobians(*(refEl.getNodes()), &refEl);
        std::vector<EMatrix> invJacs = Operator::calcInvJacobians(jacs);
        std::vector<EMatrix> diffTensors(refEl.getNumNodes(), 2.0*EMatrix::Identity(i+1, i+1));
        diffOp.setDiffTensor(diffTensors);
        std::vector<double> dV = Operator::calcMeasure(Operator::calcDetJacobians(jacs), &refEl);
        diffOp.assemble(dV, invJacs);
        LaplaceModel lapModel(&refEl);
        lapModel.allocate(1);
        lapModel.setElementNodes(refEl.getNodes());
        std::map<std::string, std::vector<double> > fm;
        fm["DiffusionTensor"] = std::vector<double>((i+1)*(i+1)*refEl.getNumNodes(), 0.0);
        for(int k = 0; k < refEl.getNumNodes(); k++){
          for(int l = 0; l < (i+1); l++){
            fm["DiffusionTensor"][(k*(i+1) + l)*(i+1) + l] = 2.0;
          }
        }
        lapModel.setFieldMap(&fm);
        lapModel.compute();
        CHECK(*(lapModel.getLocalMatrix()) == *(diffOp.getMatrix()));
        CHECK(*(lapModel.getLocalRHS()) == EVector::Zero(refEl.getNumNodes()));
      };
    }
  }

};


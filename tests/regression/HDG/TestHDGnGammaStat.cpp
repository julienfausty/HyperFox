#include <catch2/catch.hpp>
#include <string>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <boost/filesystem.hpp>
#include <fstream>
#include "Mesh.h"
#include "Field.h"
#include "ZoltanPartitioner.h"
#include "HDF5Io.h"
#include "HDGnGammaModel.h"
#include "IntegratedDirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "NonLinearWrapper.h"
#include "TestUtils.h"

using namespace hfox;

using namespace nGamma;

namespace hdg{

void manufacturedNGammaStat(std::vector<double> x, std::vector<double> kn, std::vector<double> kGam, std::vector<double> * nGamma){
  nGamma->at(0) = 2.0 + std::cos(kn[0]*x[0] + kn[1]*x[1]);
  nGamma->at(1) = std::cos(kGam[0]*x[0] + kGam[1]*x[1]);
};//manufacturedNGammaStat

void manufacturedMagneticNGammaStat(std::vector<double> x, std::vector<double> * b){
  b->at(0) = -x[1]; b->at(1) = x[0];
};//manufacturedMagneticNGammaStat

double manufacturedSourceNGammaStat(std::vector<double> x, std::vector<double> kn, std::vector<double> kGam, double diffCoeff, int i){
  double res = 0.0;
  if(i == 0){
    res = diffCoeff*(std::pow(kn[0], 2) + std::pow(kn[1], 2))*std::cos(kn[0]*x[0] + kn[1]*x[1]) - (x[0]*kGam[1] - x[1]*kGam[0])*std::sin(kGam[0]*x[0] + kGam[1]*x[1]);
  } else if (i == 1){
    res = diffCoeff*(std::pow(kGam[0], 2) + std::pow(kGam[1], 2))*std::cos(kGam[0]*x[0] + kGam[1]*x[1]) - (x[0]*kn[1] - x[1]*kn[0])*std::sin(kn[0]*x[0] + kn[1]*x[1]);
    res += -2.0*((std::sin(kGam[0]*x[0] + kGam[1]*x[1])*std::cos(kGam[0]*x[0] + kGam[1]*x[1]))/(2.0 + std::cos(kn[0]*x[0] + kn[1]*x[1])))*(x[0]*kGam[1] - x[1]*kGam[0]);
    res += ((std::pow(std::cos(kGam[0]*x[0] + kGam[1]*x[1]),2)*std::sin(kn[0]*x[0] + kn[1]*x[1]))/std::pow(2.0 + std::cos(kn[0]*x[0] + kn[1]*x[1]), 2))*(x[0]*kn[1] - x[1]*kn[0]);
  }
  return res;
};//manufacturedSourceNGammaStat

};//hdg

void runHDGnGammaStat(SimRun * thisRun){
  //setup
  //Load mesh
  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order;
  thisRun->meshLocation += meshName + ".h5";
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
  //setup partitioner
  ZoltanOpts zOpts;
  zOpts.debugLevel = "0";
  ZoltanPartitioner zPart(&myMesh, zOpts);
  zPart.initialize();
  //initialize fields
  int nNodesPerEl = myMesh.getReferenceElement()->getNumNodes();
  int nNodesPerFace = myMesh.getReferenceElement()->getFaceElement()->getNumNodes();
  int nodeDim = myMesh.getNodeSpaceDimension();
  Field dirichlet(&myMesh, Face, nNodesPerFace, nodeDim);
  Field sol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field buffSol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field flux(&myMesh, Cell, nNodesPerEl, std::pow(nodeDim, 2));
  Field trace(&myMesh, Face, nNodesPerFace, nodeDim);
  Field partition(&myMesh, Node, 1, 1);
  Field b(&myMesh, Node, 1, nodeDim);
  Field D(&myMesh, Node, 1, 1);
  Field G(&myMesh, Node, 1, 1);
  Field tau(&myMesh, Face, nNodesPerFace, std::pow(nodeDim, 2)*2);
  Field anaSol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field residual(&myMesh, Cell, nNodesPerEl, 1);
  //create fieldMap
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["BufferSolution"] = &buffSol;
  fieldMap["Flux"] = &flux;
  fieldMap["Trace"] = &trace;
  fieldMap["Tau"] = &tau;
  fieldMap["Dirichlet"] = &dirichlet;
  fieldMap["b"] = &b;
  fieldMap["D"] = &D;
  fieldMap["G"] = &G;
  fieldMap["Analytical"] = &anaSol;
  //initialize models, solvers, etc.
  IntegratedDirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGnGammaModel transportMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.maxits = 10000;
  myOpts.rtol = 1e-16;
  myOpts.verbose = false;
  myOpts.solverType = KSPGMRES;
  myOpts.preconditionnerType = PCBJACOBI;
  PetscInterface petsciface(myOpts);
  HDGSolverOpts solveOpts;
  solveOpts.type = IMPLICIT;
  solveOpts.verbosity = false;
  solveOpts.doubleValuedTau = true;
  HDGSolver mySolver;
  mySolver.setOptions(solveOpts);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&transportMod);
  mySolver.setBoundaryModel(&dirMod);
  NonLinearWrapper wrapper;
  wrapper.setVerbosity(false);
  int maxNRIter = 10;
  wrapper.setMaxIterations(maxNRIter);
  wrapper.setResidualTolerance(1e-6);
  wrapper.setSolutionFields(&sol, &buffSol);
  wrapper.setSolver(&mySolver);
  //setup outputs
  std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/nGammaStat/";
  writeDir += meshName;
  //if(zPart.getRank() == 0){
    //boost::filesystem::create_directory(writeDir);
  //}
  hdfio.setField("Solution", &sol);
  hdfio.setField("Flux", &flux);
  hdfio.setField("b", &b);
  hdfio.setField("Analytical", &anaSol);
  hdfio.setField("Residual", &residual);
  hdfio.setField("Partition", &partition);
  //define scalars
  double diffCoeff = 1.0;//diffusive coeff
  double carLen = 1.0;
  std::vector<double> kn = {3.0*M_PI, 3.0*M_PI};
  std::vector<double> kGam = {3.0*M_PI, 3.0*M_PI};
  //create buffers  
  std::vector<int> cell;
  std::vector<double> node;
  std::vector<int> face2Cell;
  std::vector<int> ibuffer;
  std::vector< std::vector<int> > iibuffer;
  std::vector<double> dbuffer(nodeDim, 0.0);
  //initialize field values
  for(int i = 0; i < myMesh.getNumberCells(); i++){
    myMesh.getCell(i, &cell);
    for(int j = 0; j < cell.size(); j++){
      int locInd = zPart.global2LocalNode(cell[j]);
      if(locInd != -1){
        myMesh.getPoint(locInd, &node);
      } else {
        myMesh.getGhostPoint(cell[j], &node);
      }
      hdg::manufacturedNGammaStat(node, kn, kGam, &dbuffer);
      for(int k = 0; k < dbuffer.size(); k++){
        anaSol.getValues()->at((i*nNodesPerEl + j)*nodeDim + k) = dbuffer[k];
        sol.getValues()->at((i*nNodesPerEl + j)*nodeDim + k) = (double)(k == 0);
      }
    }
  }
  std::copy(sol.getValues()->begin(), sol.getValues()->end(), buffSol.getValues()->begin());
  std::fill(D.getValues()->begin(), D.getValues()->end(), diffCoeff);
  std::fill(G.getValues()->begin(), G.getValues()->end(), diffCoeff);
  for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
    myMesh.getPoint(iNode, &node);
    hdg::manufacturedMagneticNGammaStat(node, &dbuffer);
    for(int k = 0; k < dbuffer.size(); k++){
      b.getValues()->at(iNode*nodeDim + k) = dbuffer[k];
    }
  }
  for(int iFace = 0; iFace < myMesh.getNumberFaces(); iFace++){
    myMesh.getFace(iFace, &cell);
    for(int iN = 0; iN < nNodesPerFace; iN++){
      for(int iCell = 0; iCell < 2; iCell++){
        EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          = EMatrix::Identity(nodeDim, nodeDim)*(diffCoeff/carLen)*2.0;
      }
    }
  }
  //compute boundary field
  for(std::set<int>::const_iterator itset = myMesh.getBoundaryFaces()->begin(); itset != myMesh.getBoundaryFaces()->end(); itset++){
    int locFace = zPart.global2LocalFace(*itset);
    if(locFace != -1){
      myMesh.getFace(locFace, &cell);
      for(int k = 0; k < cell.size(); k++){
        int locNode = zPart.global2LocalNode(cell[k]);
        if(locNode != -1){
          myMesh.getPoint(locNode, &node);
        } else {
          myMesh.getGhostPoint(cell[k], &node);
        }
        hdg::manufacturedNGammaStat(node, kn, kGam, &dbuffer);
        for(int dof = 0; dof < nodeDim; dof++){
          dirichlet.getValues()->at((locFace * nNodesPerFace + k)*nodeDim + dof) = dbuffer[dof];
        }
      }
    }
  }
  //create field list
  std::vector<Field*> fieldList;
  for(auto itMap = fieldMap.begin(); itMap != fieldMap.end(); itMap++){
    fieldList.push_back(itMap->second);
  }
  fieldList.push_back(&residual);
  fieldList.push_back(&partition);
  //partition the mesh and fields
  zPart.setFields(fieldList);
  zPart.computePartition();
  zPart.update();
  //first output
  //hdfio.write(writeDir + "/res_0.h5");
  //allocating and last set ups
  mySolver.initialize();
  mySolver.allocate();
  transportMod.setSourceFunction([kn, kGam, diffCoeff](const std::vector<double> & x, int i){return hdg::manufacturedSourceNGammaStat(x, kn, kGam, diffCoeff, i);});
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
  //solve non-linear problem
  start = std::chrono::high_resolution_clock::now();
  wrapper.solve();
  end = std::chrono::high_resolution_clock::now();
  thisRun->resolution += end - start;
  //output
  //get linalg err
  const KSP * ksp = petsciface.getKSP();
  double linAlgErr;
  KSPGetResidualNorm(*ksp, &(linAlgErr));
  //calc l2Err
  double sumRes = 0, sumAna = 0;
  for(int k = 0; k < myMesh.getNumberCells(); k++){
    for(int n = 0; n < nNodesPerEl; n++){
      residual.getValues()->at(k*nNodesPerEl + n) = 0.0;
      for(int dof = 0; dof < nodeDim; dof++){
        residual.getValues()->at(k*nNodesPerEl + n) += std::pow(anaSol.getValues()->at((k*nNodesPerEl + n)*nodeDim + dof) - sol.getValues()->at((k*nNodesPerEl + n)*nodeDim + dof), 2);
      }
      residual.getValues()->at(k*nNodesPerEl + n) = std::sqrt(residual.getValues()->at(k*nNodesPerEl + n));
      sumRes += std::pow(residual.getValues()->at(k*nNodesPerEl + n), 2);
      sumAna += std::pow(anaSol.getValues()->at(k*nNodesPerEl + n), 2);
    }
  }
  zPart.updateSharedInformation();
  double l2Err = std::sqrt(TestUtils::l2ProjectionCellField(&residual, &residual, &myMesh));
  double dL2Err = std::sqrt(sumRes/sumAna);
  thisRun->linAlgErr = linAlgErr;
  thisRun->l2Err = l2Err;
  thisRun->dL2Err = dL2Err;
  //hdfio.write(writeDir + "/res_" + std::to_string(maxNRIter) + ".h5");
  end = std::chrono::high_resolution_clock::now();
  thisRun->post += end - start;
};//runHDGnGamma

TEST_CASE("Testing regression cases for the stationary HDGnGammaModel", "[regression][HDG][nGammaStat]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  meshSizes["2"] = {"1e-1"};
  //meshSizes["2"] = {"1e-1", "7e-2", "5e-2"};
  //meshSizes["2"] = {"7e-2"};
  std::vector<std::string> orders = {"1", "2", "3", "4", "5"};
  //std::vector<std::string> orders = {"3"};
  std::vector<SimRun> simRuns;
  for(auto it = meshSizes.begin(); it != meshSizes.end(); it++){
    for(auto itMs = it->second.begin(); itMs != it->second.end(); itMs++){
      for(int o = 0; o < orders.size(); o++){
        SimRun thisRun;
        thisRun.dim = it->first;
        thisRun.meshSize = *itMs;
        thisRun.order = orders[o];
        simRuns.push_back(thisRun);
      }
    }
  }

  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/nGammaStat/";
  std::string writeFile = "Breakdown.csv";
  int nParts;
  MPI_Comm_size(MPI_COMM_WORLD, &nParts);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //std::ofstream f;
  //if(rank == 0){
    //f.open(writePath + writeFile);
    //f << "dim,order,h,linAlgErr,l2Err,dL2Err,avgRuntime,maxRuntime,minRuntime\n" << std::flush;
  //}
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGnGammaStat(&(*it));
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    CHECK(it->l2Err < 1);
    double timeBuff, maxRuntime, avgRuntime, minRuntime;
    timeBuff = it->runtime.count();
    MPI_Allreduce(&timeBuff, &maxRuntime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&timeBuff, &minRuntime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&timeBuff, &avgRuntime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    avgRuntime /= nParts;
    //if(rank == 0){
      //f << it->dim << ",";
      //f << it->order << ",";
      //f << it->meshSize << ",";
      //f << it->linAlgErr << ",";
      //f << it->l2Err << ",";
      //f << it->dL2Err << ",";
      //f << avgRuntime << ",";
      //f << maxRuntime << ",";
      //f << minRuntime << "\n";
      //f << std::flush;
    //}
  }
  //if(rank == 0){
    //f << std::endl;
    //f.close();
  //}
};

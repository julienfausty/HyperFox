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
#include "HDGConnectionLaplacianModel.h"
#include "RungeKutta.h"
#include "IntegratedDirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "IP2NodeSolver.h"
#include "TestUtils.h"

using namespace hfox;



struct ConnectionLaplacianRun{
  double R = 0.0;
  double r = 0.0;
  std::string h = "";
  std::string o = "";
  std::string writeDir = "";
  double thetaFreq = 0.0;
  double phiFreq = 0.0;
  double l2Err = 0.0;
  double linAlgErr = 0.0;
  double T = 0.0;
};

void computeTorusAngles(const std::vector<double> & xyz, std::vector<double> * phitheta, ConnectionLaplacianRun * run){
  phitheta->resize(2, 0.0);
  phitheta->at(1) = std::atan2(xyz[2]/run->r, std::sqrt(std::pow(xyz[0],2) + std::pow(xyz[1], 2)) - run->R);
  phitheta->at(0) = std::atan2(xyz[1]/(run->R + run->r*std::cos(phitheta->at(1))), xyz[0]/(run->R + run->r*std::cos(phitheta->at(1))));
};//computeTorusAngles

double computeAnaSolTorus(const std::vector<double> & xyz, ConnectionLaplacianRun * run){
  std::vector<double> phiTheta(2, 0.0);
  computeTorusAngles(xyz, &phiTheta, run);
  return std::cos(run->phiFreq*phiTheta[0] + run->thetaFreq*phiTheta[1]);
};//computeAnaSolTorus

double computeSourceTorus(const std::vector<double> & xyz, ConnectionLaplacianRun * run){
  std::vector<double> phiTheta(2, 0.0);
  computeTorusAngles(xyz, &phiTheta, run);
  double res = run->thetaFreq*std::sin(phiTheta[1])*std::sin(run->phiFreq * phiTheta[0] + run->phiFreq * phiTheta[1])/((run->R + run->r*std::cos(phiTheta[1]))*run->r);
  res -= (std::pow(run->thetaFreq/run->r, 2) + std::pow(run->phiFreq/(run->R + run->r*std::cos(phiTheta[1])), 2))*std::cos(run->phiFreq*phiTheta[0] + run->thetaFreq*phiTheta[1]);
  return -res;
};

void runConnectionLaplacianTorus(ConnectionLaplacianRun * run){
  //setup
  //Load mesh
  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  std::string inputDir = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "correctedQuarterTorus_h" + run->h + "_p" + run->o + ".h5";
  Mesh myMesh(2, std::stoi(run->o), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(inputDir + meshName);
  //setup partitioner
  ZoltanOpts zOpts;
  zOpts.debugLevel = "0";
  ZoltanPartitioner zPart(&myMesh, zOpts);
  zPart.initialize();
  //initialize fields
  const ReferenceElement * refEl = myMesh.getReferenceElement();
  int nodeDim = myMesh.getNodeSpaceDimension();
  int refDim = refEl->getDimension();
  int nNodesEl = refEl->getNumNodes();
  int nNodesFc = refEl->getFaceElement()->getNumNodes();
  Field sol(&myMesh, Cell, nNodesEl, 1);
  Field flux(&myMesh, Cell, nNodesEl, nodeDim);
  Field trace(&myMesh, Face, nNodesFc, 1);
  Field tau(&myMesh, Face, nNodesFc, 1);
  Field jacobian(&myMesh, Cell, nNodesEl, refDim*nodeDim);
  Field metric(&myMesh, Cell, nNodesEl, refDim*refDim);
  Field source(&myMesh, Cell, nNodesEl, 1);
  Field anaSol(&myMesh, Cell, nNodesEl, 1);
  Field residual(&myMesh, Cell, nNodesEl, 1);
  Field partition(&myMesh, Node, 1, 1);
  Field dirichlet(&myMesh, Face, nNodesFc, 1);
  Field dirichletFlag(&myMesh, Node, 1, 1);
  //create field map
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["Flux"] = &flux;
  fieldMap["Trace"] = &trace;
  fieldMap["Tau"] = &tau;
  fieldMap["Jacobian"] = &jacobian;
  fieldMap["Metric"] = &metric;
  fieldMap["Dirichlet"] = &dirichlet;
  //initialize models, solvers, etc.
  HDGConnectionLaplacianModel mod(refEl);
  mod.setEmbeddingDimension(3);
  IntegratedDirichletModel bMod(myMesh.getReferenceElement()->getFaceElement());
  std::set<int> dfaces;
  PetscOpts myOpts;
  myOpts.maxits = 3000;
  myOpts.rtol = 1e-12;
  myOpts.verbose = false;
  myOpts.preconditionnerType = PCLU;
  PetscInterface petsciface(myOpts);
  HDGSolverOpts solveOpts;
  solveOpts.verbosity = false;
  HDGSolver mySolver;
  mySolver.setOptions(solveOpts);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&mod);
  //mySolver.setBoundaryCondition(&bMod, &dfaces);
  mySolver.setBoundaryModel(&bMod);
  //setup outputs
  std::string output = "quarterTorus_h" + run->h + "_p" + run->o + "_sol.h5";
  hdfio.setField("Solution", &sol);
  hdfio.setField("Flux", &flux);
  hdfio.setField("Jacobian", &jacobian);
  hdfio.setField("Metric", &metric);
  hdfio.setField("Source", &source);
  hdfio.setField("Dirichlet", &dirichletFlag);
  hdfio.setField("Analytical", &anaSol);
  hdfio.setField("Residual", &residual);
  hdfio.setField("Partition", &partition);
  //setup fields
  std::vector<int> cell(nNodesEl, 0);
  std::vector<double> node(3, 0.0);
  std::fill(tau.getValues()->begin(), tau.getValues()->end(), 1.0);
  for(int iEl = 0; iEl < myMesh.getNumberCells(); iEl++){
    myMesh.getCell(iEl, &cell);
    for(int iN = 0; iN < nNodesEl; iN++){
      myMesh.getPoint(cell[iN], &node);
      anaSol.getValues()->at(iEl*nNodesEl + iN) = computeAnaSolTorus(node, run);
      source.getValues()->at(iEl*nNodesEl + iN) = computeSourceTorus(node, run);
    }
  }
  IP2NodeSolver ip2n(&myMesh);
  ip2n.setField(&jacobian);
  ip2n.setFunction(IP2NodeSolver::jacobianVals);
  ip2n.solve();
  ip2n.setField(&metric);
  ip2n.setFunction(IP2NodeSolver::metricVals);
  ip2n.solve();
  //Boundary conditions
  //double tol = std::stod(run->h)/10.0;
  double tol = 1e-4;
  std::set<int>::const_iterator it;
  for(it = myMesh.getBoundaryFaces()->begin(); it != myMesh.getBoundaryFaces()->end(); it++){
    myMesh.getFace(*it, &cell);
    for(int iN = 0; iN < nNodesFc; iN++){
      myMesh.getPoint(cell[iN], &node);
      dirichlet.getValues()->at((*it)*nNodesFc + iN) = computeAnaSolTorus(node, run);
    }
  }
  //create field list
  std::vector<Field*> fieldList;
  for(auto itMap = fieldMap.begin(); itMap != fieldMap.end(); itMap++){
    fieldList.push_back(itMap->second);
  }
  fieldList.push_back(&anaSol);
  fieldList.push_back(&residual);
  fieldList.push_back(&partition);
  fieldList.push_back(&source);
  fieldList.push_back(&dirichletFlag);
  //partition the mesh and fields
  zPart.setFields(fieldList);
  zPart.computePartition();
  zPart.update();
  //setup partition field
  int rank = zPart.getRank();
  for(int iN = 0; iN < myMesh.getNumberPoints(); iN++){
    partition.getValues()->at(iN) = rank;
  }
  //allocate
  mySolver.initialize();
  mySolver.allocate();
  mod.setSourceFunction([run](const std::vector<double> & node){return computeSourceTorus(node, run);});
  //solve
  mySolver.assemble();
  mySolver.solve();
  //compute errors
  //get linalg err
  const KSP * ksp = petsciface.getKSP();
  double linAlgErr;
  KSPGetResidualNorm(*ksp, &(linAlgErr));
  //calc l2Err
  for(int k = 0; k < myMesh.getNumberCells(); k++){
    for(int n = 0; n < nNodesEl; n++){
      residual.getValues()->at(k*nNodesEl + n) = anaSol.getValues()->at(k*nNodesEl + n) - sol.getValues()->at(k*nNodesEl + n);
    }
  }
  double l2Err = std::sqrt(TestUtils::l2ProjectionCellField(&residual, &residual, &myMesh));
  run->linAlgErr = linAlgErr;
  run->l2Err = l2Err;
  zPart.updateSharedInformation();
  //output
  hdfio.write(run->writeDir + output);
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  run->T = (end - start).count();
};//runConnectionLaplacianTorus

TEST_CASE("Testing regression cases for the HDGConnectionLaplacianModel", "[regression][HDG][ConnectionLaplacianModel]"){
  //Large radius
  double R = 3.0;
  //Small radius
  double r = 1.0;
  //ThetaFreq
  double nfreq = 3.0;
  //PhiFreq
  double mfreq = 3.0;
  //Mesh sizes
  std::vector<std::string> hs = {"1e-0", "7e-1", "5e-1", "2e-1", "1e-1"};
  //Polynomial orders
  std::vector<std::string> orders = {"1", "2", "3", "4", "5"};
  //Output
  std::string writeDir = "/home/julien/workspace/RenaissanceFusion/Postprocess/results/TorusRegression/";
  //Runs
  std::vector<ConnectionLaplacianRun> runs;
  for(int iH = 0; iH < hs.size(); iH++){
    for(int iO = 0; iO < orders.size(); iO++){
      ConnectionLaplacianRun run;
      run.R = R;
      run.r = r;
      run.h = hs[iH];
      run.o = orders[iO];
      run.thetaFreq = nfreq;
      run.phiFreq = mfreq;
      run.writeDir = writeDir;
      runs.push_back(run);
    }
  }
  for(int iR = 0; iR < runs.size(); iR++){
    std::cout << "Running h = " << runs[iR].h << ", p = " << runs[iR].o << std::endl;
    runConnectionLaplacianTorus(&(runs[iR]));
  }
  int nParts;
  MPI_Comm_size(MPI_COMM_WORLD, &nParts);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ofstream f;
  if(rank == 0){
    f.open(writeDir + "Breakdown.csv");
    f << "R,r,h,o,l2Err,linAlgErr,runtime" << std::flush;
    for(auto it = runs.begin(); it < runs.end(); it++){
      f << "\n";
      f << it->R << ",";
      f << it->r << ",";
      f << it->h << ",";
      f << it->o << ",";
      f << it->l2Err << ",";
      f << it->linAlgErr << ",";
      f << it->T << std::flush;
    }
    f.close();
  }
};

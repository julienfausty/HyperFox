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
#include "Operator.h"
#include "HDGBurgersModel.h"
#include "HDGDiffusionSource.h"
#include "RungeKutta.h"
#include "IntegratedDirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "NonLinearWrapper.h"
#include "TestUtils.h"

using namespace hfox;

namespace hdg{


void analyticalBurgersStatx5(const std::vector<double> & x, std::vector<double> * sol){
  sol->resize(x.size(), 0.0);
  double buff;
  sol->at(0) = std::pow(x[0], 5);
  //sol->at(0) = 1.0;
};

double sourceBurgersStat(const std::vector<double> & x, int i, double D){
  double res = 0.0;
  double cBuff, dBuff;
  double ccBuff;
  if(i == 0){
    //res = 2.0 * std::pow(x[0], 1);
    res = 5.0 * std::pow(x[0], 9) - D * 5.0 * 4.0 * std::pow(x[0], 3);
    //res = - D * 5 * 4 * std::pow(x[0], 3);
  }
  return res;
};

double potentialiBurgersStat(std::vector<double> x, double & A, std::vector<double> & x0){
  return (0.001*A*std::exp((1+x0[0])*A)*(1+x[0]) + (std::exp(A*(x[0] - x0[0])) + std::exp(-A*(x[0]-x0[0])))*std::cos(A*(x[1] - x0[1])));
};

void potentialiGradBurgersStat(std::vector<double> x, double & A, std::vector<double> & x0, std::vector<double> * grad){
  grad->resize(x.size());
  grad->at(0) =  0.001*A*std::exp((1+x0[0])*A) + A * (std::exp(A*(x[0] - x0[0])) - std::exp(-A*(x[0]-x0[0])))*std::cos(A*(x[1] - x0[1]));
  grad->at(1) = -A * (std::exp(A*(x[0] - x0[0])) + std::exp(-A*(x[0]-x0[0])))*std::sin(A*(x[1] - x0[1]));
};

void potentialiHessBurgersStat(std::vector<double> x, double & A, std::vector<double> & x0, std::vector<double> * hess){
  hess->resize(std::pow(x.size(), 2));
  std::vector<double> grad;
  potentialiGradBurgersStat(x, A, x0, &grad);
  double pot = potentialiBurgersStat(x, A, x0);
  hess->at(0) = std::pow(A, 2.0) * pot;
  hess->at(1) = -A * A * (std::exp(A*(x[0] - x0[0])) - std::exp(-A*(x[0]-x0[0])))*std::sin(A*(x[1] - x0[1]));
  hess->at(2) = hess->at(1);
  hess->at(3) = -std::pow(A, 2)*pot;
};

double potentialBurgersStat(std::vector<double> x, double & A, std::vector<double> & x0){
  return (potentialiBurgersStat(x, A, x0));
};

void potentialGradBurgersStat(std::vector<double> x, double & A, std::vector<double> & x0, std::vector<double> * grad){
  grad->resize(x.size());
  potentialiGradBurgersStat(x, A, x0, grad);
};

void potentialHessBurgersStat(std::vector<double> x, double & A, std::vector<double> & x0, std::vector<double> * hess){
  potentialiHessBurgersStat(x, A, x0, hess);
};

void analyticalBurgersStat(std::vector<double> x, double D, std::vector<double> * sol){
  double A = 0.5;
  std::vector<double> x0 = {1.0, 0.0};
  sol->resize(x.size(), 0.0);
  //sol->at(0) = 1.0;
  potentialGradBurgersStat(x, A, x0, sol);
  double buff = potentialBurgersStat(x, A, x0);
  EMap<EVector>(sol->data(), sol->size()) *= -2.0*D/buff;
};

void analyticalBurgersGradStat(std::vector<double> x, double D, std::vector<double> * gradSol){
  double A = 0.5;
  std::vector<double> x0 = {1.0, 0.0};
  gradSol->resize(std::pow(x.size(), 2), 0.0);
  double pot = potentialBurgersStat(x, A, x0);
  std::vector<double> gradPot;
  potentialGradBurgersStat(x, A, x0, &gradPot);
  potentialHessBurgersStat(x, A, x0, gradSol);
  EMap<EMatrix> res(gradSol->data(), x.size(), x.size());
  res *= -2.0*D/pot;
  EMap<EVector> gPot(gradPot.data(), gradPot.size());
  res += (2.0*D/std::pow(pot, 2.0))*(gPot * gPot.transpose());
};

};//hdg

void runHDGBurgersStat(SimRun * thisRun, HDGSolverType globType){
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
  Field diffCoeff(&myMesh, Node, 1, 1);
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
  fieldMap["DiffusionTensor"] = &diffCoeff;
  fieldMap["Analytical"] = &anaSol;
  //initialize models, solvers, etc.
  IntegratedDirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGBurgersModel transportMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.maxits = 10000;
  myOpts.rtol = 1e-16;
  myOpts.verbose = false;
  myOpts.solverType = KSPGMRES;
  myOpts.preconditionnerType = PCBJACOBI;
  PetscInterface petsciface(myOpts);
  HDGSolverOpts solveOpts;
  solveOpts.type = globType;
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
  wrapper.setVerbosity(true);
  int maxNRIter = 10;
  wrapper.setVerbosity(false);
  wrapper.setMaxIterations(maxNRIter);
  wrapper.setResidualTolerance(1e-6);
  wrapper.setSolutionFields(&sol, &buffSol);
  wrapper.setSolver(&mySolver);
  //setup outputs
  //std::string writeDir = "/home/jfausty/workspace/Postprocess/results/BurgersStat/HDG/";
  std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/BurgersStat/HDG/";
  writeDir += meshName;
  //if(zPart.getRank() == 0){
    //boost::filesystem::create_directory(writeDir);
  //}
  hdfio.setField("Solution", &sol);
  hdfio.setField("Buffer", &buffSol);
  hdfio.setField("Flux", &flux);
  hdfio.setField("Analytical", &anaSol);
  hdfio.setField("Residual", &residual);
  hdfio.setField("Partition", &partition);
  //define scalars
  double D = 1.0;//diffusive coeff
  double carLen = std::stod(thisRun->meshSize);
  int dimGrad = std::pow(nodeDim, 2);
  //create buffers  
  std::vector<int> cell;
  std::vector<double> node;
  std::vector<int> face2Cell;
  std::vector<int> ibuffer;
  std::vector< std::vector<int> > iibuffer;
  std::vector<double> dbuffer(nodeDim, 0.0);
  std::vector<EVector> normals(nNodesPerFace, EVector::Zero(nodeDim));
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
      hdg::analyticalBurgersStat(node, D, &dbuffer);
      //hdg::analyticalBurgersStatx5(node, &dbuffer);
      for(int k = 0; k < dbuffer.size(); k++){
        anaSol.getValues()->at((i*nNodesPerEl + j)*nodeDim + k) = dbuffer[k];
        //sol.getValues()->at((i*nNodesPerEl + j)*nodeDim + k) = dbuffer[k];
      }
      //sol.getValues()->at((i*nNodesPerEl + j)*nodeDim) *= 0.5;
      //hdg::analyticalBurgersGradStat(node, D, &dbuffer);
      //for(int k = 0; k < dbuffer.size(); k++){
        //flux.getValues()->at((i*nNodesPerEl + j)*dimGrad + k) = dbuffer[k];
      //}
    }
  }
  std::fill(diffCoeff.getValues()->begin(), diffCoeff.getValues()->end(), D);
  for(int iFace = 0; iFace < myMesh.getNumberFaces(); iFace++){
    //normals = TestUtils::calculateOutwardNormal(&myMesh, iFace);
    myMesh.getFace(iFace, &cell);
    for(int j = 0; j < cell.size(); j++){
      int locInd = zPart.global2LocalNode(cell[j]);
      if(locInd != -1){
        myMesh.getPoint(locInd, &node);
      } else {
        myMesh.getGhostPoint(cell[j], &node);
      }
      hdg::analyticalBurgersStat(node, D, &dbuffer);
      //hdg::analyticalBurgersStatx5(node, &dbuffer);
      for(int k = 0; k < dbuffer.size(); k++){
        //trace.getValues()->at((iFace*nNodesPerFace + j)*nodeDim + k) = dbuffer[k];
      }
    }
    //myMesh.getFace2Cell(iFace, &face2Cell);
    //anaSol.getSliceValues(face2Cell, &dbuffer);
    //myMesh.getSliceCells(face2Cell, &iibuffer);
    for(int iN = 0; iN < nNodesPerFace; iN++){
      for(int iCell = 0; iCell < 2; iCell++){
        //int offset = iCell*nNodesPerEl + std::distance(iibuffer[iCell].begin(), std::find(iibuffer[iCell].begin(), iibuffer[iCell].end(), cell[iN]));
        //double pun = normals[iN].dot(EMap<EVector>(dbuffer.data() + offset*nodeDim, nodeDim));
        //EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          //= EMatrix::Identity(nodeDim, nodeDim)*(std::fabs(pun) + D);
        EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          = EMatrix::Identity(nodeDim, nodeDim)*(D/carLen)*2.0;
      }
    }
  }
  //create field list
  std::vector<Field*> fieldList;
  for(auto itMap = fieldMap.begin(); itMap != fieldMap.end(); itMap++){
    fieldList.push_back(itMap->second);
  }
  //fieldList.push_back(&anaSol);
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
  //transportMod.setSourceFunction([D](const std::vector<double> & x, int i){return hdg::sourceBurgersStat(x, i, D);});
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
  //compute necessary fields
  std::copy(sol.getValues()->begin(), sol.getValues()->end(), buffSol.getValues()->begin());
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
        hdg::analyticalBurgersStat(node, D, &dbuffer);
        //hdg::analyticalBurgersStatx5(node, &dbuffer);
        for(int dof = 0; dof < nodeDim; dof++){
          dirichlet.getValues()->at((locFace * nNodesPerFace + k)*nodeDim + dof) = dbuffer[dof];
        }
      }
    }
  }

  zPart.updateSharedInformation();
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
  //double l2Err = TestUtils::l2AnalyticalErrCellField(&sol, [D](std::vector<double> & x, std::vector<double> * sol){hdg::analyticalBurgersStat(x, D, sol);}, &myMesh);
  double dL2Err = std::sqrt(sumRes/sumAna);
  thisRun->linAlgErr = linAlgErr;
  thisRun->l2Err = l2Err;
  thisRun->dL2Err = dL2Err;
  //hdfio.write(writeDir + "/res_" + std::to_string(maxNRIter) + ".h5");
  end = std::chrono::high_resolution_clock::now();
  thisRun->post += end - start;
};



TEST_CASE("Testing stationary regression cases for the HDGBurgersModel", "[regression][HDG][BurgersModelStat]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  meshSizes["2"] = {"2e-1", "1e-1"};
  //meshSizes["2"] = {"1e-1"};
  std::vector<std::string> orders = {"1", "2", "3"};
  //std::vector<std::string> orders = {"1"};
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

  //std::string writePath = "/home/jfausty/workspace/Postprocess/results/BurgersStat/HDG/";
  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/BurgersStat/HDG/";
  HDGSolverType globType = IMPLICIT;
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
    runHDGBurgersStat(&(*it), globType);
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

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
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "NonLinearWrapper.h"
#include "TestUtils.h"

using namespace hfox;

namespace hdg{

double potentialiBurgersStat(std::vector<double> x, double & Ax, double & Ay, std::vector<double> & x0){
  return ((std::exp(Ax*(x[0] - x0[0])) + std::exp(-Ax*(x[0]-x0[0])))*std::cos(Ay*(x[1] - x0[1])));
};

void potentialiGradBurgersStat(std::vector<double> x, double & Ax, double & Ay, std::vector<double> & x0, std::vector<double> * grad){
  grad->resize(x.size());
  grad->at(0) = Ax * (std::exp(Ax*(x[0] - x0[0])) - std::exp(-Ax*(x[0]-x0[0])))*std::cos(Ay*(x[1] - x0[1]));
  grad->at(1) = -Ay * (std::exp(Ax*(x[0] - x0[0])) + std::exp(-Ax*(x[0]-x0[0])))*std::sin(Ay*(x[1] - x0[1]));
};

void potentialiHessBurgersStat(std::vector<double> x, double & Ax, double & Ay, std::vector<double> & x0, std::vector<double> * hess){
  hess->resize(std::pow(x.size(), 2));
  std::vector<double> grad;
  potentialiGradBurgersStat(x, Ax, Ay, x0, &grad);
  double pot = potentialiBurgersStat(x, Ax, Ay, x0);
  hess->at(0) = std::pow(Ax, 2.0) * pot;
  hess->at(1) = -Ay * Ax * (std::exp(Ax*(x[0] - x0[0])) - std::exp(-Ax*(x[0]-x0[0])))*std::sin(Ay*(x[1] - x0[1]));
  hess->at(2) = hess->at(1);
  hess->at(3) = -std::pow(Ay, 2)*pot;
};

double potentialBurgersStat(std::vector<double> x, double & Ax0, double & Ay0, double & Ax1, double & Ay1, std::vector<double> & x0){
  //return (potentialiBurgersStat(x, Ax0, Ay0, x0) + potentialiBurgersStat(x, Ax1, Ay1, x0));
  return (potentialiBurgersStat(x, Ax0, Ay0, x0));
};

void potentialGradBurgersStat(std::vector<double> x, double & Ax0, double & Ay0, double & Ax1, double & Ay1, std::vector<double> & x0, std::vector<double> * grad){
  grad->resize(x.size());
  potentialiGradBurgersStat(x, Ax0, Ay0, x0, grad);
  //std::vector<double> buff;
  //potentialiGradBurgersStat(x, Ax1, Ay1, x0, &buff);
  //EMap<EVector>(grad->data(), grad->size()) += EMap<EVector>(buff.data(), buff.size());
};

void potentialHessBurgersStat(std::vector<double> x, double & Ax0, double & Ay0, double & Ax1, double & Ay1, std::vector<double> & x0, std::vector<double> * hess){
  potentialiHessBurgersStat(x, Ax0, Ay0, x0, hess);
  //std::vector<double> buff;
  //potentialiHessBurgersStat(x, Ax1, Ay1, x0, &buff);
  //EMap<EMatrix>(hess->data(), x.size(), x.size()) += EMap<EMatrix>(buff.data(), x.size(), x.size());
};

void analyticalBurgersStat(std::vector<double> x, double D, std::vector<double> * sol){
  double Ax0 = 0.1;
  double Ay0 = 0.1;
  double Ax1 = 2.0;
  double Ay1 = 1.0;
  std::vector<double> x0 = {0.5, 0.5};
  sol->resize(x.size(), 0.0);
  //sol->at(0) = 1.0;
  potentialGradBurgersStat(x, Ax0, Ay0, Ax1, Ay1, x0, sol);
  double buff = potentialBurgersStat(x, Ax0, Ay0, Ax1, Ay1, x0);
  if(buff != 0){
    EMap<EVector>(sol->data(), sol->size()) *= -2.0*D/buff;
  } else {
    EMap<EVector>(sol->data(), sol->size()) *= 0.0;
  }
};

void analyticalBurgersGradStat(std::vector<double> x, double D, std::vector<double> * gradSol){
  double Ax0 = 0.1;
  double Ay0 = 0.1;
  double Ax1 = 2.0;
  double Ay1 = 1.0;
  std::vector<double> x0 = {0.5, 0.5};
  gradSol->resize(std::pow(x.size(), 2), 0.0);
  double pot = potentialBurgersStat(x, Ax0, Ay0, Ax1, Ay1, x0);
  std::vector<double> gradPot;
  potentialGradBurgersStat(x, Ax0, Ay0, Ax1, Ay1, x0, &gradPot);
  potentialHessBurgersStat(x, Ax0, Ay0, Ax1, Ay1, x0, gradSol);
  EMap<EMatrix> res(gradSol->data(), x.size(), x.size());
  if(pot != 0){
    res *= -2.0*D/pot;
    EMap<EVector> gPot(gradPot.data(), gradPot.size());
    res += (2.0*D/std::pow(pot, 2.0))*(gPot * gPot.transpose());
  } else {
    res = EMatrix::Zero(x.size(), x.size());
  }
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
  //initialize models, solvers, etc.
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGBurgersModel transportMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.maxits = 3000;
  myOpts.rtol = 1e-12;
  myOpts.verbose = false;
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
  int maxNRIter = 50;
  wrapper.setMaxIterations(maxNRIter);
  wrapper.setSolutionFields(&sol, &buffSol);
  wrapper.setSolver(&mySolver);
  //setup outputs
  //std::string writeDir = "/home/jfausty/workspace/Postprocess/results/BurgersStat/HDG/";
  std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/BurgersStat/HDG/";
  writeDir += meshName;
  if(zPart.getRank() == 0){
    boost::filesystem::create_directory(writeDir);
  }
  hdfio.setField("Solution", &sol);
  hdfio.setField("Buffer", &buffSol);
  hdfio.setField("Flux", &flux);
  hdfio.setField("Analytical", &anaSol);
  hdfio.setField("Residual", &residual);
  hdfio.setField("Partition", &partition);
  //define scalars
  double D = 1.0;//1e-2;//diffusive coeff
  double carLen = 1.0;//std::sqrt(D*timeStep);
  int dimGrad = std::pow(nodeDim, 2);
  //create buffers  
  std::vector<int> cell;
  std::vector<double> node;
  std::vector<int> face2Cell;
  std::vector<int> ibuffer;
  std::vector<double> dbuffer;
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
    for(int iN = 0; iN < nNodesPerFace; iN++){
      for(int iCell = 0; iCell < 2; iCell++){
        EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          = EMatrix::Identity(nodeDim, nodeDim);
      }
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
  //partition the mesh and fields
  zPart.setFields(fieldList);
  zPart.computePartition();
  zPart.update();
  //first output
  hdfio.write(writeDir + "/res_0.h5");
  //allocating and last set ups
  mySolver.initialize();
  mySolver.allocate();
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
  //mySolver.assemble();
  //mySolver.solve();
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
  thisRun->linAlgErr += linAlgErr;
  thisRun->l2Err += l2Err;
  thisRun->dL2Err += dL2Err;
  hdfio.write(writeDir + "/res_" + std::to_string(maxNRIter) + ".h5");
  end = std::chrono::high_resolution_clock::now();
  thisRun->post += end - start;
};



TEST_CASE("Testing stationary regression cases for the HDGBurgersModel", "[regression][HDG][BurgersModelStat]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  meshSizes["2"] = {"1e-1", "7e-2", "5e-2", "2e-2"};
  //meshSizes["2"] = {"2e-2"};
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

  //std::string writePath = "/home/jfausty/workspace/Postprocess/results/BurgersStat/HDG/";
  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/BurgersStat/HDG/";
  HDGSolverType globType = IMPLICIT;
  std::string writeFile = "Breakdown.csv";
  int nParts;
  MPI_Comm_size(MPI_COMM_WORLD, &nParts);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ofstream f;
  if(rank == 0){
    f.open(writePath + writeFile);
    f << "dim,order,h,linAlgErr,l2Err,dL2Err,avgRuntime,maxRuntime,minRuntime\n" << std::flush;
  }
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
    if(rank == 0){
      f << it->dim << ",";
      f << it->order << ",";
      f << it->meshSize << ",";
      f << it->linAlgErr << ",";
      f << it->l2Err << ",";
      f << it->dL2Err << ",";
      f << avgRuntime << ",";
      f << maxRuntime << ",";
      f << minRuntime << "\n";
      f << std::flush;
    }
  }
  if(rank == 0){
    f << std::endl;
    f.close();
  }
};

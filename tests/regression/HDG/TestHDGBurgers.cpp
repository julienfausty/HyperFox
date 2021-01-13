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

void manufacturedBurgersSolution(double t, std::vector<double> x, double D, std::vector<double> * sol){
  sol->resize(x.size(), 0.0);
  std::fill(sol->begin(), sol->end(), sol->size());
  sol->at(0) = std::pow(x[0], 5) + std::pow(t, 3);
  sol->at(1) = std::pow(x[1], 4) + std::pow(t, 2);
}

void manufacturedBurgersSolutionGrad(double t, std::vector<double> x, double D, std::vector<double> * sol){
  sol->resize(std::pow(x.size(), 2), 0.0);
  std::fill(sol->begin(), sol->end(), sol->size());
  sol->at(0) = 5.0 * std::pow(x[0], 4);
  sol->at(3) = 4.0 * std::pow(x[1], 3);
}

double manufacturedSourceBurgers(const double t, const std::vector<double> & x, int i, double D){
  double res = 0.0;
  if(i == 0){
    res = (std::pow(x[0], 5) + std::pow(t, 3))*5.0 * std::pow(x[0], 4) - D * 5.0 * 4.0 * std::pow(x[0], 3) + 3.0*std::pow(t, 2.0);
  } else if(i == 1){
    res = (std::pow(x[1], 4) + std::pow(t, 2))*4.0 * std::pow(x[1], 3) - D * 4.0 * 3.0 * std::pow(x[1], 2) + 2.0*std::pow(t, 1.0);
  }
  return res;
};

double potentialiBurgers(double t, std::vector<double> x, double & Ax, double & Ay, std::vector<double> & x0){
  return (std::exp((std::pow(Ax, 2) - std::pow(Ay, 2))*t)*(std::exp(Ax*(x[0] - x0[0])) + std::exp(-Ax*(x[0]-x0[0])))*std::cos(Ay*(x[1] - x0[1])));
};

void potentialiGradBurgers(double t, std::vector<double> x, double & Ax, double & Ay, std::vector<double> & x0, std::vector<double> * grad){
  grad->resize(x.size());
  grad->at(0) = Ax * std::exp((std::pow(Ax, 2) - std::pow(Ay, 2))*t)*(std::exp(Ax*(x[0] - x0[0])) - std::exp(-Ax*(x[0]-x0[0])))*std::cos(Ay*(x[1] - x0[1]));
  grad->at(1) = - Ay * std::exp((std::pow(Ax, 2) - std::pow(Ay, 2))*t)*(std::exp(Ax*(x[0] - x0[0])) + std::exp(-Ax*(x[0]-x0[0])))*std::sin(Ay*(x[1] - x0[1]));
};

void potentialiHessBurgers(double t, std::vector<double> x, double & Ax, double & Ay, std::vector<double> & x0, std::vector<double> * hess){
  hess->resize(std::pow(x.size(), 2));
  std::vector<double> grad;
  potentialiGradBurgers(t, x, Ax, Ay, x0, &grad);
  double pot = potentialiBurgers(t, x, Ax, Ay, x0);
  hess->at(0) = std::pow(Ax, 2.0) * pot;
  hess->at(1) = -Ay * Ax * std::exp((std::pow(Ax, 2) - std::pow(Ay, 2))*t)*(std::exp(Ax*(x[0] - x0[0])) - std::exp(-Ax*(x[0]-x0[0])))*std::sin(Ay*(x[1] - x0[1]));
  hess->at(2) = hess->at(1);
  hess->at(3) = -std::pow(Ay, 2)*pot;
};

double potentialBurgers(double t, std::vector<double> x, double & Ax0, double & Ay0, double & Ax1, double & Ay1, std::vector<double> & x0, std::vector<double> & x1){
  return (potentialiBurgers(t, x, Ax0, Ay0, x0) + potentialiBurgers(t, x, Ax1, Ay1, x1));
  //return (potentialiBurgers(t, x, Ax0, Ay0, x0));
};

void potentialGradBurgers(double t, std::vector<double> x, double & Ax0, double & Ay0, double & Ax1, double & Ay1, std::vector<double> & x0, std::vector<double> & x1, std::vector<double> * grad){
  grad->resize(x.size());
  potentialiGradBurgers(t, x, Ax0, Ay0, x0, grad);
  std::vector<double> buff;
  potentialiGradBurgers(t, x, Ax1, Ay1, x1, &buff);
  EMap<EVector>(grad->data(), grad->size()) += EMap<EVector>(buff.data(), buff.size());
};

void potentialHessBurgers(double t, std::vector<double> x, double & Ax0, double & Ay0, double & Ax1, double & Ay1, std::vector<double> & x0, std::vector<double> & x1, std::vector<double> * hess){
  potentialiHessBurgers(t, x, Ax0, Ay0, x0, hess);
  std::vector<double> buff;
  potentialiHessBurgers(t, x, Ax1, Ay1, x1, &buff);
  EMap<EMatrix>(hess->data(), x.size(), x.size()) += EMap<EMatrix>(buff.data(), x.size(), x.size());
};

void analyticalBurgers(double t, std::vector<double> x, double D, std::vector<double> * sol){
  double Ax0 = 1.2;
  double Ay0 = 1.8;
  double Ax1 = 0.6;
  double Ay1 = 1.2;
  std::vector<double> x0 = {0.25, 0.25};
  std::vector<double> x1 = {0.75, 0.75};
  sol->resize(x.size());
  //sol->at(0) = 1.0;
  potentialGradBurgers(t, x, Ax0, Ay0, Ax1, Ay1, x0, x1, sol);
  double buff = potentialBurgers(t, x, Ax0, Ay0, Ax1, Ay1, x0, x1);
  if(buff != 0){
    EMap<EVector>(sol->data(), sol->size()) *= -2.0*D/buff;
  } else {
    EMap<EVector>(sol->data(), sol->size()) *= 0.0;
  }
};

void analyticalBurgersGrad(double t, std::vector<double> x, double D, std::vector<double> * gradSol){
  double Ax0 = 1.2;
  double Ay0 = 1.8;
  double Ax1 = 0.6;
  double Ay1 = 1.2;
  std::vector<double> x0 = {0.25, 0.25};
  std::vector<double> x1 = {0.75, 0.75};
  gradSol->resize(std::pow(x.size(), 2), 0.0);
  double pot = potentialBurgers(t, x, Ax0, Ay0, Ax1, Ay1, x0, x1);
  std::vector<double> gradPot;
  potentialGradBurgers(t, x, Ax0, Ay0, Ax1, Ay1, x0, x1, &gradPot);
  potentialHessBurgers(t, x, Ax0, Ay0, Ax1, Ay1, x0, x1, gradSol);
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

void runHDGBurgers(SimRun * thisRun, HDGSolverType globType){
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
  //setup time scheme
  RKType rkType;
  std::string rkStr = thisRun->rk;
  if(rkStr == "BEuler"){
    rkType = BEuler;
  }else if(rkStr == "ALX2"){
    rkType = ALX2;
  }else if(rkStr == "IMidpoint"){
    rkType = IMidpoint;
  }else if(rkStr == "RK43"){
    rkType = RK43;
  }else if(rkStr == "FEuler"){
    rkType = FEuler;
  }else if(rkStr == "SSPRK3"){
    rkType = SSPRK3;
  } else{
    rkType = BEuler;
  }
  std::vector<std::string> auxiliaries = {"Flux", "Trace"};
  RungeKutta ts(myMesh.getReferenceElement(), rkType, auxiliaries);
  int nStages = ts.getNumStages();
  //initialize fields
  int nNodesPerEl = myMesh.getReferenceElement()->getNumNodes();
  int nNodesPerFace = myMesh.getReferenceElement()->getFaceElement()->getNumNodes();
  int nodeDim = myMesh.getNodeSpaceDimension();
  Field dirichlet(&myMesh, Face, nNodesPerFace, nodeDim);
  Field sol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field buffSol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field oldSol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field flux(&myMesh, Cell, nNodesPerEl, std::pow(nodeDim, 2));
  Field oldFlux(&myMesh, Cell, nNodesPerEl, std::pow(nodeDim, 2));
  Field trace(&myMesh, Face, nNodesPerFace, nodeDim);
  Field oldTrace(&myMesh, Face, nNodesPerFace, nodeDim);
  Field partition(&myMesh, Node, 1, 1);
  Field diffCoeff(&myMesh, Node, 1, 1);
  Field tau(&myMesh, Face, nNodesPerFace, std::pow(nodeDim, 2)*2);
  Field anaSol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field residual(&myMesh, Cell, nNodesPerEl, 1);
  std::vector<Field> rkStages(nStages, Field(&myMesh, Cell, nNodesPerEl, nodeDim));
  std::vector<Field> rkFluxStages(nStages, Field(&myMesh, Cell, nNodesPerEl, std::pow(nodeDim, 2)));
  std::vector<Field> rkTraceStages(nStages, Field(&myMesh, Face, nNodesPerFace, nodeDim));
  //create fieldMap
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["BufferSolution"] = &buffSol;
  fieldMap["OldSolution"] = &oldSol;
  fieldMap["Flux"] = &flux;
  fieldMap["OldFlux"] = &oldFlux;
  fieldMap["Trace"] = &trace;
  fieldMap["OldTrace"] = &oldTrace;
  for(int i = 0; i < nStages; i++){
    fieldMap["RKStage_" + std::to_string(i)] = &(rkStages[i]);
    fieldMap["RKStage_Flux_" + std::to_string(i)] = &(rkFluxStages[i]);
    fieldMap["RKStage_Trace_" + std::to_string(i)] = &(rkTraceStages[i]);
  }
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
  wrapper.setVerbosity(false);
  wrapper.setMaxIterations(12);
  wrapper.setSolutionFields(&sol, &buffSol);
  wrapper.setSolver(&mySolver);
  //setup outputs
  std::string writeDir = "/home/jfausty/workspace/Postprocess/results/Burgers/HDG/";
  //std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/Burgers/HDG/";
  if(globType == WEXPLICIT){
    writeDir += "WExp/";
  } else if(globType == SEXPLICIT){
    writeDir += "SExp/";
  } else {
    writeDir += "Imp/";
  }
  if(rkStr == "BEuler"){
    writeDir += "BEuler/";
  }else if(rkStr == "ALX2"){
    writeDir += "ALX2/";
  }else if(rkStr == "IMidpoint"){
    writeDir += "IMidpoint/";
  }else if(rkStr == "RK43"){
    writeDir += "RK43/";
  }else if(rkStr == "FEuler"){
    writeDir += "FEuler/";
  }else if(rkStr == "SSPRK3"){
    writeDir += "SSPRK3/";
  } else{
    writeDir += "Misc/";
  }
  writeDir += meshName + "_dt-" + thisRun->timeStep;
  //if(zPart.getRank() == 0){
    //boost::filesystem::create_directory(writeDir);
  //}
  hdfio.setField("Solution", &sol);
  hdfio.setField("Flux", &flux);
  hdfio.setField("Analytical", &anaSol);
  hdfio.setField("Residual", &residual);
  hdfio.setField("Partition", &partition);
  //define scalars
  double D = 1e-1;//diffusive coeff
  double timeStep = std::stod(thisRun->timeStep);
  double carLen = std::sqrt(D*timeStep);
  double t = 0;
  int dimGrad = std::pow(nodeDim, 2);
  double timeEnd = 1.0;
  int nIters = timeEnd/timeStep;
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
      //hdg::analyticalBurgers(t, node, D, &dbuffer);
      hdg::manufacturedBurgersSolution(t, node, D, &dbuffer);
      for(int k = 0; k < dbuffer.size(); k++){
        sol.getValues()->at((i*nNodesPerEl + j)*nodeDim + k) = dbuffer[k];
      }
      //hdg::analyticalBurgersGrad(t, node, D, &dbuffer);
      hdg::manufacturedBurgersSolutionGrad(t, node, D, &dbuffer);
      for(int k = 0; k < dbuffer.size(); k++){
        flux.getValues()->at((i*nNodesPerEl + j)*dimGrad + k) = dbuffer[k];
      }
    }
  }
  std::copy(sol.getValues()->begin(), sol.getValues()->end(), anaSol.getValues()->begin());
  std::fill(diffCoeff.getValues()->begin(), diffCoeff.getValues()->end(), D);
  for(int iFace = 0; iFace < myMesh.getNumberFaces(); iFace++){
    for(int iN = 0; iN < nNodesPerFace; iN++){
      for(int iCell = 0; iCell < 2; iCell++){
        EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          = EMatrix::Identity(nodeDim, nodeDim);
        //EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          //= EMatrix::Identity(nodeDim, nodeDim)*(D/carLen);
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
  //hdfio.write(writeDir + "/res_0.h5");
  //allocating and last set ups
  ts.setTimeStep(timeStep);
  transportMod.setTimeScheme(&ts);
  mySolver.initialize();
  mySolver.allocate();
  transportMod.setSourceFunction([&t, D](const std::vector<double> & x, int i){return hdg::manufacturedSourceBurgers(t, x, i, D);});
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
  //time iteration
  int i = 0;
  //ProgressBar pbar;
  //pbar.setIterIndex(&i);
  //pbar.setNumIterations(nIters);
  //std::cout << "Simulation (d=" + thisRun->dim + ", h=" + thisRun->meshSize + ", p=" + thisRun->order + ", dt=" + thisRun->timeStep + ", rank=" + std::to_string(zPart.getRank()) + ")" << std::endl;
  //pbar.update();
  for(i = 0; i < nIters; i++){
    t += timeStep;
    //compute necessary fields
    for(int j = 0; j < myMesh.getNumberCells(); j++){
      myMesh.getCell(j, &cell);
      for(int l = 0; l < cell.size(); l++){
        int locInd = zPart.global2LocalNode(cell[l]);
        if(locInd != -1){
          myMesh.getPoint(locInd, &node);
        } else {
          myMesh.getGhostPoint(cell[l], &node);
        }
        //hdg::analyticalBurgers(t, node, D, &dbuffer);
        hdg::manufacturedBurgersSolution(t, node, D, &dbuffer);
        for(int k = 0; k < dbuffer.size(); k++){
          anaSol.getValues()->at((j*nNodesPerEl + l)*nodeDim + k) = dbuffer[k];
        }
      }
    }
    std::copy(sol.getValues()->begin(), sol.getValues()->end(), buffSol.getValues()->begin());
    std::copy(sol.getValues()->begin(), sol.getValues()->end(), oldSol.getValues()->begin());
    std::copy(flux.getValues()->begin(), flux.getValues()->end(), oldFlux.getValues()->begin());
    std::copy(trace.getValues()->begin(), trace.getValues()->end(), oldTrace.getValues()->begin());
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
          //hdg::analyticalBurgers(t, node, D, &dbuffer);
          hdg::manufacturedBurgersSolution(t, node, D, &dbuffer);
          for(int dof = 0; dof < nodeDim; dof++){
            dirichlet.getValues()->at((locFace * nNodesPerFace + k)*nodeDim + dof) = dbuffer[dof];
          }
        }
      }
    }
    //need to try upwind tau
    //for(int iFace = 0; iFace < myMesh.getNumberFaces(); iFace++){
      //normals = TestUtils::calculateOutwardNormal(&myMesh, iFace);
      //myMesh.getFace(iFace, &cell);
      //myMesh.getFace2Cell(iFace, &face2Cell);
      //for(int iCell = 0; iCell < face2Cell.size(); iCell++){
        //int locInd = zPart.global2LocalElement(face2Cell[iCell]);
        //if(locInd != -1){
          //sol.getValues(locInd, &dbuffer);
          //myMesh.getCell(locInd, &ibuffer);
          //for(int j = 0; j < nNodesPerFace; j++){
            //int offset = std::distance(ibuffer.begin(), std::find(ibuffer.begin(), ibuffer.end(), cell[j]));
            //double pun = normals[j].dot(EMap<EVector>(dbuffer.data() + offset*nodeDim, nodeDim));
            ////double mult = 0.0;
            //EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + j)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim) 
              //= EMatrix::Identity(nodeDim, nodeDim);
            ////EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + j)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim) 
            ////= EMatrix::Identity(nodeDim, nodeDim)*mult;
            ////EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + j)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim) 
            ////= EMatrix::Identity(nodeDim, nodeDim)*(std::fabs(pun) + D/carLen);
            ////if((pun > 0 and iCell == 0) or (pun < 0 and iCell == 1)){
            ////EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + j)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim) 
            ////= EMatrix::Identity(nodeDim, nodeDim)*(std::fabs(pun) + D/carLen);
            //////EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + j)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim) 
            //////= EMatrix::Identity(nodeDim, nodeDim)*(std::fabs(pun));
            ////} else {
            ////EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + j)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim) 
            ////= EMatrix::Identity(nodeDim, nodeDim)*(D/carLen);
            //////EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + j)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim) 
            //////= EMatrix::Identity(nodeDim, nodeDim)*0.0;
            ////}
          //}
        //}
      //}
    //}
    zPart.updateSharedInformation();
    //solve non-linear problem
    start = std::chrono::high_resolution_clock::now();
    for(int k = 0; k < ts.getNumStages(); k++){
      wrapper.solve();
      //mySolver.assemble();
      //mySolver.solve();
      ts.computeStage(&fieldMap);
    }
    ts.computeSolution(&fieldMap);
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
    //std::cout << "l2Err at time " << t << ": " << l2Err << std::endl;
    if(i != (nIters-1)){
      thisRun->linAlgErr += linAlgErr*timeStep;
      thisRun->l2Err += l2Err*timeStep;
      thisRun->dL2Err += dL2Err*timeStep;
    } else{
      thisRun->linAlgErr += linAlgErr*timeStep/2;
      thisRun->l2Err += l2Err*timeStep/2;
      thisRun->dL2Err += dL2Err*timeStep/2;
    }
    double quot = t/(5e-3);
    //double quot = 0.0;
    double rem = quot - ((int)quot);
    //std::cout << "rem: " << rem << std::endl;
    //if(rem < timeStep/(5e-3)){
      //hdfio.write(writeDir + "/res_" + std::to_string(i+1) + ".h5");
    //}
    end = std::chrono::high_resolution_clock::now();
    thisRun->post += end - start;
    //pbar.update();
    if(thisRun->l2Err > 1.0){
      break;
    }
  }
};



TEST_CASE("Testing regression cases for the HDGBurgersModel", "[regression][HDG][BurgersModel]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  meshSizes["2"] = {"2e-1", "1e-1"};
  //meshSizes["2"] = {"7e-2"};
  //std::vector<std::string> timeSteps = {"1e-2", "5e-3", "2e-3", "1e-3", "5e-4", "2e-4", "1e-4", "5e-5", "2e-5"};
  std::vector<std::string> timeSteps = {"1e-2"};
  //std::vector<std::string> timeSteps = {"1e-3"};
  std::vector<std::string> orders = {"1", "2", "3"};
  //std::vector<std::string> orders = {"3"};
  std::vector<std::string> rkTypes = {"BEuler"};
  std::vector<SimRun> simRuns;
  for(auto it = meshSizes.begin(); it != meshSizes.end(); it++){
    for(auto itMs = it->second.begin(); itMs != it->second.end(); itMs++){
      for(int o = 0; o < orders.size(); o++){
        for(int dt = 0; dt < timeSteps.size(); dt++){
          for(int rki = 0; rki < rkTypes.size(); rki++){
            SimRun thisRun;
            thisRun.dim = it->first;
            thisRun.meshSize = *itMs;
            thisRun.order = orders[o];
            thisRun.timeStep = timeSteps[dt];
            thisRun.rk = rkTypes[rki];
            simRuns.push_back(thisRun);
          }
        }
      }
    }
  }

  std::string writePath = "/home/jfausty/workspace/Postprocess/results/Burgers/HDG/";
  //std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/Burgers/HDG/";
  //HDGSolverType globType = WEXPLICIT;
  HDGSolverType globType = IMPLICIT;
  //HDGSolverType globType = SEXPLICIT;
  std::string writeFile = "Breakdown.csv";
  if(globType == WEXPLICIT){
    writePath += "WExp/";
  } else if(globType == SEXPLICIT){
    writePath += "SExp/";
  } else {
    writePath += "Imp/";
  }
  if(rkTypes[0] == "BEuler"){
    writePath += "BEuler/";
  }else if(rkTypes[0] == "ALX2"){
    writePath += "ALX2/";
  }else if(rkTypes[0] == "IMidpoint"){
    writePath += "IMidpoint/";
  }else if(rkTypes[0] == "RK43"){
    writePath += "RK43/";
  }else if(rkTypes[0] == "FEuler"){
    writePath += "FEuler/";
  }else if(rkTypes[0] == "SSPRK3"){
    writePath += "SSPRK3/";
  } else{
    writePath += "Misc/";
  }
  int nParts;
  MPI_Comm_size(MPI_COMM_WORLD, &nParts);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //std::ofstream f;
  //if(rank == 0){
    //f.open(writePath + writeFile);
    //f << "dim,order,h,dt,linAlgErr,l2Err,dL2Err,avgRuntime,maxRuntime,minRuntime\n" << std::flush;
  //}
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGBurgers(&(*it), globType);
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
      //f << it->timeStep << ",";
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

#include <catch2/catch.hpp>
#include <string>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <boost/filesystem.hpp>
#include <fstream>
#include "Mesh.h"
#include "ZoltanPartitioner.h"
#include "Field.h"
#include "HDF5Io.h"
#include "Operator.h"
#include "HDGConvectionDiffusionReactionSource.h"
#include "RungeKutta.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "TestUtils.h"

using namespace hfox;

double analyticalCDRSHDG(const double t, const std::vector<double> & x, double & v, std::vector<double> & cvel, double & D){
  double sigma = 0.1;
  std::vector<double> c = {0.3, 0.3};
  double vt = v*t;
  double twoSigSq = 2.0*std::pow(sigma, 2);
  double fourKapT = 4.0*D*t;
  double xbar = (x[0] - cvel[0])*std::cos(vt) + (x[1] - cvel[1])*std::sin(vt);
  double ybar = -(x[0] - cvel[0])*std::sin(vt) + (x[1] - cvel[1])*std::cos(vt);
  double res = (twoSigSq/(twoSigSq + fourKapT))*std::exp(-(std::pow(xbar - c[0] + cvel[0], 2) + std::pow(ybar - c[1] + cvel[1], 2))/(twoSigSq + fourKapT));
  return res;
};

double analyticalCDRSHDGGrad(const double t, const std::vector<double> & x, double & v, std::vector<double> & cvel, double & D, int k){
  double sigma = 0.1;
  std::vector<double> c = {0.3, 0.3};
  double vt = v*t;
  double twoSigSq = 2.0*std::pow(sigma, 2);
  double fourKapT = 4.0*D*t;
  double cvt = std::cos(vt);
  double svt = std::sin(vt);
  double xbar = (x[0] - cvel[0])*cvt + (x[1] - cvel[1])*svt;
  double ybar = -(x[0] - cvel[0])*svt + (x[1] - cvel[1])*cvt;
  double res = -(2.0*twoSigSq/std::pow(twoSigSq + fourKapT, 2))*std::exp(-(std::pow(xbar - c[0] + cvel[0], 2) + std::pow(ybar - c[1] + cvel[1], 2))/(twoSigSq + fourKapT));
  if(k == 0){
    res *= (cvt*(xbar - c[0] + cvel[0]) - svt*(ybar - c[1] + cvel[1]));
  } else if(k == 1){
    res *= (svt*(xbar - c[0] + cvel[0]) + cvt*(ybar - c[1] + cvel[1]));
  } else {
    res *= 0;
  }
  return res;
};

void runHDGCDRS(SimRun * thisRun,  HDGSolverType globType){
  //setup
  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order;
  thisRun->meshLocation += meshName + ".h5";
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
  ZoltanOpts zOpts;
  zOpts.debugLevel = "0";
  ZoltanPartitioner zPart(&myMesh, zOpts);
  zPart.initialize();
  //std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/parallel/CDRS/";
  std::string writeDir = "/home/jfausty/workspace/Postprocess/results/parallel/CDRS/";
  if(globType == WEXPLICIT){
    writeDir += "WExp/";
  } else if(globType == SEXPLICIT){
    writeDir += "SExp/";
  } else {
    writeDir += "Imp/";
  }
  RKType rkType;
  std::string rkStr = thisRun->rk;
  if(rkStr == "BEuler"){
    writeDir += "BEuler/";
    rkType = BEuler;
  }else if(rkStr == "ALX2"){
    writeDir += "ALX2/";
    rkType = ALX2;
  }else if(rkStr == "IMidpoint"){
    writeDir += "IMidpoint/";
    rkType = IMidpoint;
  }else if(rkStr == "RK43"){
    writeDir += "RK43/";
    rkType = RK43;
  }else if(rkStr == "FEuler"){
    writeDir += "FEuler/";
    rkType = FEuler;
  }else if(rkStr == "SSPRK3"){
    writeDir += "SSPRK3/";
    rkType = SSPRK3;
  } else{
    writeDir += "Misc/";
    rkType = BEuler;
  }
  writeDir += std::to_string(zPart.getNumPartitions()) + "/";
  writeDir += meshName + "_dt-" + thisRun->timeStep;
  //if(zPart.getRank() == 0){
    //boost::filesystem::create_directory(writeDir);
  //}
  std::vector<std::string> auxiliaries = {"Flux", "Trace"};
  RungeKutta ts(myMesh.getReferenceElement(), rkType, auxiliaries);
  int nStages = ts.getNumStages();
  int nNodes = myMesh.getReferenceElement()->getNumNodes();
  int nNodesPerFace = myMesh.getReferenceElement()->getFaceElement()->getNumNodes();
  int nodeDim = myMesh.getNodeSpaceDimension();
  Field dirichlet(&myMesh, Face, nNodesPerFace, 1);
  Field sol(&myMesh, Cell, nNodes, 1);
  Field oldSol(&myMesh, Cell, nNodes, 1);
  std::vector<Field> rkStages(nStages, Field(&myMesh, Cell, nNodes, 1));
  Field flux(&myMesh, Cell, nNodes, nodeDim);
  Field oldFlux(&myMesh, Cell, nNodes, nodeDim);
  std::vector<Field> rkFluxStages(nStages, Field(&myMesh, Cell, nNodes, nodeDim));
  Field trace(&myMesh, Face, nNodesPerFace, 1);
  Field oldTrace(&myMesh, Face, nNodesPerFace, 1);
  std::vector<Field> rkTraceStages(nStages, Field(&myMesh, Face, nNodesPerFace, 1));
  Field vel(&myMesh, Node, 1, nodeDim);
  Field partition(&myMesh, Node, 1, 1);
  std::vector<double> node(nodeDim);
  double timeStep = std::stod(thisRun->timeStep);
  double v = 4.0; std::vector<double> cvel = {0.5, 0.5};
  double D = 1e-2; 
  double carLen = std::sqrt(D*timeStep);
  for(int i = 0; i < myMesh.getNumberPoints(); i++){
    myMesh.getPoint(i, &node);
    vel.getValues()->at(i*nodeDim) = -v*(node[1] - cvel[1]);
    vel.getValues()->at(i*nodeDim + 1) = v*(node[0] - cvel[0]);
  }
  Field diff(&myMesh, Node, 1, 1);
  std::fill(diff.getValues()->begin(), diff.getValues()->end(), D);
  Field tau(&myMesh, Face, nNodesPerFace, 2);
  tau.setDoubleValued(true);
  Field anaSol(&myMesh, Cell, nNodes, 1);
  Field residual(&myMesh, Cell, nNodes, 1);
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
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
  fieldMap["Velocity"] = &vel;
  fieldMap["DiffusionTensor"] = &diff;
  hdfio.setField("Solution", &sol);
  hdfio.setField("Partition", &partition);
  const std::set<int> * boundary = myMesh.getBoundaryFaces();
  std::vector<int> cell(nNodesPerFace, 0);
  std::vector< std::vector<double> > nodes(nNodesPerFace, std::vector<double>(nodeDim, 0.0));
  double t = 0;
  std::vector<double> locV(nodeDim, 0.0);
  std::vector<EVector> normals(nNodesPerFace, EVector::Zero(nodeDim));
  double projV = 0.0;
  for(int i = 0; i < myMesh.getNumberFaces(); i++){
    normals = TestUtils::calculateOutwardNormal(&myMesh, i);
    myMesh.getFace(i, &cell);
    for(int j = 0; j < nNodesPerFace; j++){
      vel.getValues(cell[j], &locV);
      projV = normals[j].dot(EMap<EVector>(locV.data(), locV.size()));
      //tau.getValues()->at(2*(i*nNodesPerFace + j)) = 0.0;
      //tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = 0.0;
      tau.getValues()->at(2*(i*nNodesPerFace + j)) = std::fabs(projV) + D/carLen;
      tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = std::fabs(projV) + D/carLen;
      //if(projV > 0){
        //tau.getValues()->at(2*(i*nNodesPerFace + j)) = std::fabs(projV) + D/carLen;
        ////tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = 0.0; //Imp
        //tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = D/carLen; //Exp
      //} else {
        ////tau.getValues()->at(2*(i*nNodesPerFace + j)) = 0.0; //Imp
        //tau.getValues()->at(2*(i*nNodesPerFace + j)) = D/carLen; //Exp
        //tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = std::fabs(projV) + D/carLen;
      //}
    }
  }
  for(int i = 0; i < myMesh.getNumberCells(); i++){
    myMesh.getCell(i, &cell);
    for(int j = 0; j < nNodes; j++){
      myMesh.getPoint(cell[j], &node);
      sol.getValues()->at(i*nNodes + j) = analyticalCDRSHDG(t, node, v, cvel, D);
      for(int k = 0; k < nodeDim; k++){
        flux.getValues()->at((i*nNodes + j)*nodeDim + k) = analyticalCDRSHDGGrad(t, node, v, cvel, D, k);
      }
    }
  }
  for(int i = 0; i < myMesh.getNumberFaces(); i++){
    myMesh.getFace(i, &cell);
    for(int j = 0; j < nNodesPerFace; j++){
      myMesh.getPoint(cell[j], &node);
      trace.getValues()->at(i*nNodesPerFace + j) = analyticalCDRSHDG(t, node, v, cvel, D);;
    }
  }
  std::vector<Field*> fieldList;
  for(auto itMap = fieldMap.begin(); itMap != fieldMap.end(); itMap++){
    fieldList.push_back(itMap->second);
  }
  fieldList.push_back(&anaSol);
  fieldList.push_back(&residual);
  fieldList.push_back(&partition);
  zPart.setFields(fieldList);
  zPart.computePartition();
  zPart.update();
  for(int i = 0; i < myMesh.getNumberPoints(); i++){
    partition.getValues()->at(i) = zPart.getRank();
  }
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGConvectionDiffusionReactionSource transportMod(myMesh.getReferenceElement());
  ts.setTimeStep(timeStep);
  transportMod.setTimeScheme(&ts);
  PetscOpts myOpts;
  myOpts.maxits = 20000;
  myOpts.rtol = 1e-12;
  myOpts.verbose = false;
  PetscInterface petsciface(myOpts);
  HDGSolverOpts solveOpts;
  solveOpts.type = globType;
  solveOpts.verbosity = false;
  HDGSolver mySolver;
  mySolver.setOptions(solveOpts);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&transportMod);
  mySolver.setBoundaryModel(&dirMod);
  mySolver.initialize();
  mySolver.allocate();
  //hdfio.write(writeDir + "/res_0.h5");
  double timeEnd = 0.2;
  int nIters = timeEnd / timeStep;
  //int nIters = 2;
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
  int i = 0;
  //ProgressBar pbar;
  //pbar.setIterIndex(&i);
  //pbar.setNumIterations(nIters);
  //if(zPart.getRank() == 0){
    //std::cout << "Simulation (d=" + thisRun->dim + ", h=" + thisRun->meshSize + ", p=" + thisRun->order + ", dt=" + thisRun->timeStep + ")" << std::endl;
  //}
  //pbar.update();
  for(i = 0; i < nIters; i++){
    t += timeStep;
    //analytical sol
    for(int iCell = 0; iCell < myMesh.getNumberCells(); iCell++){
      myMesh.getCell(iCell, &cell);
      for(int k = 0; k < cell.size(); k++){
        int iNode = zPart.global2LocalNode(cell[k]);
        if(iNode != -1){
          myMesh.getPoint(iNode, &node);
        } else {
          myMesh.getGhostPoint(cell[k], &node);
        }
        anaSol.getValues()->at(iCell*nNodes + k) = analyticalCDRSHDG(t, node, v, cvel, D);
      }
    }
    for(auto it = boundary->begin(); it != boundary->end(); it++){
      int locFc = zPart.global2LocalFace(*it);
      myMesh.getFace(locFc, &cell);
      for(int k = 0; k < cell.size(); k++){
        int iNode = zPart.global2LocalNode(cell[k]);
        if(iNode != -1){
          myMesh.getPoint(iNode, &node);
        } else {
          myMesh.getGhostPoint(cell[k], &node);
        }
        dirichlet.getValues()->at(locFc*nNodesPerFace+k) = analyticalCDRSHDG(t, node, v, cvel, D);
      }
    }
    //copy sol into oldsol
    for(int iSol = 0; iSol < sol.getLength(); iSol++){
      oldSol.getValues()->at(iSol) = sol.getValues()->at(iSol);
    }
    for(int iSol = 0; iSol < flux.getLength(); iSol++){
      oldFlux.getValues()->at(iSol) = flux.getValues()->at(iSol);
    }
    for(int iSol = 0; iSol < trace.getLength(); iSol++){
      oldTrace.getValues()->at(iSol) = trace.getValues()->at(iSol);
    }
    zPart.updateSharedInformation();
    for(int k = 0; k < ts.getNumStages(); k++){
      start = std::chrono::high_resolution_clock::now();
      mySolver.assemble();
      end = std::chrono::high_resolution_clock::now();
      thisRun->assembly += end - start;
      start = std::chrono::high_resolution_clock::now();
      mySolver.solve();
      ts.computeStage(&fieldMap);
      zPart.updateSharedInformation();
      end = std::chrono::high_resolution_clock::now();
      thisRun->resolution += end - start;
    }
    start = std::chrono::high_resolution_clock::now();
    ts.computeSolution(&fieldMap);
    zPart.updateSharedInformation();
    end = std::chrono::high_resolution_clock::now();
    thisRun->resolution += end - start;
    start = std::chrono::high_resolution_clock::now();
    //get linalg err
    const KSP * ksp = petsciface.getKSP();
    double linAlgErr;
    KSPGetResidualNorm(*ksp, &(linAlgErr));
    //calc l2Err
    cell.resize(nNodes);
    double sumRes = 0, sumAna = 0;
    for(int k = 0; k < myMesh.getNumberCells(); k++){
      for(int j = 0; j < nNodes; j++){
        residual.getValues()->at(k*nNodes + j) = anaSol.getValues()->at(k*nNodes + j) - sol.getValues()->at(k*nNodes + j);
        sumRes += std::pow(residual.getValues()->at(k*nNodes + j), 2);
        sumAna += std::pow(anaSol.getValues()->at(k*nNodes + j), 2);
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
    if(rem < timeStep/(5e-3)){
      //hdfio.write(writeDir + "/res_" + std::to_string(i+1) + ".h5");
    }
    end = std::chrono::high_resolution_clock::now();
    thisRun->post += end - start;
    //pbar.update();
    if(thisRun->l2Err > 1.0){
      break;
    }
  }
};

TEST_CASE("Testing regression cases for ConvectionDiffusionReactionSource", "[parallel][HDG][ConvectionDiffusionReactionSource]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  //meshSizes["2"] = {"5e-2"};
  meshSizes["2"] = {"2e-1", "1e-1"};
  //std::vector<std::string> timeSteps = {"1e-2", "5e-3", "2e-3", "1e-3", "5e-4", "2e-4", "1e-4", "5e-5", "2e-5"};
  std::vector<std::string> timeSteps = {"1e-2"};
  std::vector<std::string> orders = {"1", "2"};
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

  //std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/parallel/CDRS/";
  std::string writePath = "/home/jfausty/workspace/Postprocess/results/parallel/CDRS/";
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
  writePath += std::to_string(nParts) + "/";
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //std::ofstream f;
  //if(rank == 0){
    //f.open(writePath + writeFile);
    //f << "dim, order, h, linAlgErr, l2Err, dL2Err, avgRuntime, maxRuntime, minRuntime\n" << std::flush;
  //}
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGCDRS(&(*it), globType);
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
      //f << it->dim << ", ";
      //f << it->order << ", ";
      //f << it->meshSize << ", ";
      //f << it->linAlgErr << ", ";
      //f << it->l2Err << ", ";
      //f << it->dL2Err << ", ";
      //f << avgRuntime << ", ";
      //f << maxRuntime << ", ";
      //f << minRuntime << "\n";
      //f << std::flush;
    //}
  }
  //if(rank == 0){
    //f << std::endl;
    //f.close();
  //}
};

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
#include "RungeKutta.h"
#include "IntegratedDirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "NonLinearWrapper.h"
#include "TestUtils.h"

using namespace hfox;

using namespace nGamma;

namespace hdg{

void manufacturedNGamma(double t, std::vector<double> x, double c, std::vector<double> kn, std::vector<double> kGam, std::vector<double> * nGamma){
  nGamma->resize(2, 0.0);
  nGamma->at(0) = 2.0 + std::cos(t*c + kn[0]*x[0] + kn[1]*x[1]);
  nGamma->at(1) = std::cos(t*c + kGam[0]*x[0] + kGam[1]*x[1]);
};//manufacturedNGamma

void manufacturedMagneticNGamma(std::vector<double> x, std::vector<double> * b){
  b->resize(2, 0);
  b->at(0) = -x[1]; b->at(1) = x[0];
};//manufacturedMagneticNGamma

double manufacturedSourceNGamma(double t, std::vector<double> x, double c, std::vector<double> kn, std::vector<double> kGam, double diffCoeff, int i){
  double res = 0.0;
  if(i == 0){
    res = -c*std::sin(t*c + kn[0]*x[0] + kn[1]*x[1]) + diffCoeff*(std::pow(kn[0], 2) + std::pow(kn[1], 2))*std::cos(t*c + kn[0]*x[0] + kn[1]*x[1]) - (x[0]*kGam[1] - x[1]*kGam[0])*std::sin(t*c + kGam[0]*x[0] + kGam[1]*x[1]);
  } else if (i == 1){
    res = -c*std::sin(t*c + kn[0]*x[0] + kn[1]*x[1]) + diffCoeff*(std::pow(kGam[0], 2) + std::pow(kGam[1], 2))*std::cos(t*c + kGam[0]*x[0] + kGam[1]*x[1]) - (x[0]*kn[1] - x[1]*kn[0])*std::sin(t*c + kn[0]*x[0] + kn[1]*x[1]);
    res += -2.0*((std::sin(t*c + kGam[0]*x[0] + kGam[1]*x[1])*std::cos(t*c + kGam[0]*x[0] + kGam[1]*x[1]))/(2.0 + std::cos(t*c + kn[0]*x[0] + kn[1]*x[1])))*(x[0]*kGam[1] - x[1]*kGam[0]);
    res += ((std::pow(std::cos(t*c + kGam[0]*x[0] + kGam[1]*x[1]),2)*std::sin(t*c + kn[0]*x[0] + kn[1]*x[1]))/std::pow(2.0 + std::cos(t*c + kn[0]*x[0] + kn[1]*x[1]), 2))*(x[0]*kn[1] - x[1]*kn[0]);
  }
  return res;
};//manufacturedSourceNGamma

void nGammaSolver(RungeKutta * ts, std::map<std::string, Field*> * fieldMap, Solver * solver){
  for(int k = 0; k < ts->getNumStages(); k++){
    solver->assemble();
    solver->solve();
    ts->computeStage(fieldMap);
  }
  ts->computeSolution(fieldMap);
}//nGammaSolver

};//hdg

void runHDGnGamma(SimRun * thisRun, HDGSolverType globType){
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
  Field b(&myMesh, Node, 1, nodeDim);
  Field D(&myMesh, Node, 1, 1);
  Field G(&myMesh, Node, 1, 1);
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
  fieldMap["b"] = &b;
  fieldMap["D"] = &D;
  fieldMap["G"] = &G;
  //initialize models, solvers, etc.
  IntegratedDirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGnGammaModel transportMod(myMesh.getReferenceElement());
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
  int maxNRIter = 10;
  wrapper.setMaxIterations(maxNRIter);
  wrapper.setResidualTolerance(1e-6);
  wrapper.setSolutionFields(&sol, &buffSol);
  wrapper.setSolver(&mySolver);
  wrapper.setLinearizedSolver([&ts, &fieldMap](Solver * solver){hdg::nGammaSolver(&ts, &fieldMap, solver);});
  //setup outputs
  std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/nGamma/";
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
  if(zPart.getRank() == 0){
    boost::filesystem::create_directory(writeDir);
  }
  hdfio.setField("Solution", &sol);
  hdfio.setField("Flux", &flux);
  hdfio.setField("b", &b);
  hdfio.setField("Analytical", &anaSol);
  hdfio.setField("Residual", &residual);
  hdfio.setField("Partition", &partition);
  //define scalars
  double diffCoeff = 1e-3;//diffusive coeff
  double timeStep = std::stod(thisRun->timeStep);
  //double carLen = 1.0;//Imp
  double carLen = std::sqrt(diffCoeff*timeStep);//Exp
  std::vector<double> kn = {M_PI, M_PI};
  std::vector<double> kGam = {M_PI, M_PI};
  double c = M_PI;
  double t = 0;
  double timeEnd = 1.0;
  int nIters = timeEnd/timeStep;
  //create buffers
  std::vector<int> cell;
  std::vector<double> node;
  std::vector<int> face2Cell;
  std::vector<int> ibuffer;
  std::vector<double> dbuffer;
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
      hdg::manufacturedNGamma(t, node, c, kn, kGam, &dbuffer);
      for(int k = 0; k < dbuffer.size(); k++){
        anaSol.getValues()->at((i*nNodesPerEl + j)*nodeDim + k) = dbuffer[k];
        sol.getValues()->at((i*nNodesPerEl + j)*nodeDim + k) = dbuffer[k];
      }
    }
  }
  std::copy(sol.getValues()->begin(), sol.getValues()->end(), buffSol.getValues()->begin());
  std::fill(D.getValues()->begin(), D.getValues()->end(), diffCoeff);
  std::fill(G.getValues()->begin(), G.getValues()->end(), diffCoeff);
  for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
    myMesh.getPoint(iNode, &node);
    hdg::manufacturedMagneticNGamma(node, &dbuffer);
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
  ts.setTimeStep(timeStep);
  transportMod.setTimeScheme(&ts);
  mySolver.initialize();
  mySolver.allocate();
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
 //time iteration
  int i = 0;
  ProgressBar pbar;
  pbar.setIterIndex(&i);
  pbar.setNumIterations(nIters);
  std::cout << "Simulation (d=" + thisRun->dim + ", h=" + thisRun->meshSize + ", p=" + thisRun->order + ", dt=" + thisRun->timeStep + ", rank=" + std::to_string(zPart.getRank()) + ")" << std::endl;
  pbar.update();
  for(i = 0; i < nIters; i++){
    t += timeStep;
    transportMod.setSourceFunction([t, c, kn, kGam, diffCoeff](const std::vector<double> & x, int k){return hdg::manufacturedSourceNGamma(t, x, c, kn, kGam, diffCoeff, k);});
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
        hdg::manufacturedNGamma(t, node, c, kn, kGam, &dbuffer);
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
          hdg::manufacturedNGamma(t, node, c, kn, kGam, &dbuffer);
          for(int dof = 0; dof < nodeDim; dof++){
            dirichlet.getValues()->at((locFace * nNodesPerFace + k)*nodeDim + dof) = dbuffer[dof];
          }
        }
      }
    }
    zPart.updateSharedInformation();
    //solve non-linear problem
    start = std::chrono::high_resolution_clock::now();
    //for(int k = 0; k < ts.getNumStages(); k++){
      //wrapper.solve();//Imp
      ////mySolver.assemble();//Exp
      ////mySolver.solve();//Exp
      //ts.computeStage(&fieldMap);
    //}
    //ts.computeSolution(&fieldMap);
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
    double rem = quot - ((int)quot);
    if(rem < timeStep/(5e-3)){
      hdfio.write(writeDir + "/res_" + std::to_string(i+1) + ".h5");
    }
    end = std::chrono::high_resolution_clock::now();
    thisRun->post += end - start;
    pbar.update();
    if((thisRun->l2Err > 1.0) or (std::isnan(thisRun->l2Err))){
      break;
    }
  }

};//runHDGnGamma

TEST_CASE("Testing regression cases for the HDGnGammaModel", "[regression][HDG][nGammaModel]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  meshSizes["2"] = {"2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["2"] = {"1e-1", "7e-2", "5e-2"};
  //meshSizes["2"] = {"1e-1"};
  //std::vector<std::string> timeSteps = {"1e-2", "5e-3", "2e-3", "1e-3", "5e-4", "2e-4", "1e-4", "5e-5", "2e-5"};
  std::vector<std::string> timeSteps = {"1e-2", "5e-3", "1e-3", "5e-4", "2e-4"};
  //std::vector<std::string> timeSteps = {"1e-3"};
  std::vector<std::string> orders = {"1", "2", "3"};
  //std::vector<std::string> orders = {"3"};
  std::vector<std::string> rkTypes = {"ALX2"};
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

  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/nGamma/";
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
  std::ofstream f;
  if(rank == 0){
    f.open(writePath + writeFile);
    f << "dim,order,h,dt,linAlgErr,l2Err,dL2Err,avgRuntime,maxRuntime,minRuntime\n" << std::flush;
  }
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGnGamma(&(*it), globType);
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
      f << it->timeStep << ",";
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

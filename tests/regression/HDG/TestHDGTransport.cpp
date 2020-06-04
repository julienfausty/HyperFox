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
#include "HDF5Io.h"
#include "Operator.h"
#include "HDGTransport.h"
#include "RungeKutta.h"
#include "Euler.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "TestUtils.h"

using namespace hfox;

double analyticalTransportHDG(const double t, const std::vector<double> & x, const std::vector<double> & v){
  std::vector<double> c(x.size(), 0.3);
  double dev = 1.0/16.0;
  double freq = 8.0*M_PI;
  std::vector<double> p(x.size(), 0.0);
  EMap<EVector>(p.data(), p.size()) = EMap<const EVector>(x.data(), x.size()) - t*EMap<const EVector>(v.data(), v.size());
  return TestUtils::morelet(p, c, freq, dev);
};

double analyticalTransportHDGGrad(const double t, const std::vector<double> & x, const std::vector<double> & v, int k){
  std::vector<double> c(x.size(), 0.3);
  double dev = 1.0/16.0;
  double freq = 8.0*M_PI;
  std::vector<double> p(x.size(), 0.0);
  EMap<EVector>(p.data(), p.size()) = EMap<const EVector>(x.data(), x.size()) - t*EMap<const EVector>(v.data(), v.size());
  return TestUtils::moreletGrad(p, c, freq, dev, k);
};

std::vector<EVector> calculateOutwardNormal(Mesh * myMesh, int faceInd){
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int dimSpace = myMesh->getNodeSpaceDimension();
  int nNodes = fEl->getNumNodes();
  std::vector<int> faceIndexes;
  std::vector< std::vector<double> > faceNodes;
  const std::vector< std::vector< std::vector<double> > > * derivShapes = fEl->getDerivShapeFunctions();
  myMesh->getFace(faceInd, &faceIndexes);
  myMesh->getSlicePoints(faceIndexes, &faceNodes);
  EVector testVec = EMap<EVector>(faceNodes[0].data(), faceNodes[0].size());
  std::vector<double> node(dimSpace, 0.0);
  std::vector<int> face2Cell(2, 0);
  std::vector<int> cell(refEl->getNumNodes(), 0);
  myMesh->getFace2Cell(faceInd, &face2Cell);
  myMesh->getCell(face2Cell[0], &cell);
  std::vector<int>::iterator it;
  int innerNode = 0;
  for(int i = 0; i < cell.size(); i++){
    it = std::find(faceIndexes.begin(), faceIndexes.end(), cell[i]);
    if(it == faceIndexes.end()){
      innerNode = i;
      break;
    }
  }
  myMesh->getPoint(cell[innerNode], &node);
  testVec -= EMap<EVector>(node.data(), node.size());
  std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
  EMatrix jacobian(dimSpace, dimSpace - 1);
  for(int i = 0; i < nNodes; i++){
    jacobian *= 0.0;
    for(int j = 0; j < nNodes; j++){
      jacobian += EMap<const EVector>(faceNodes[j].data(), faceNodes[j].size())*(EMap<const EVector>(derivShapes->at(i)[j].data(), derivShapes->at(i)[j].size())).transpose();
    }
    res[i] = (jacobian*jacobian.transpose()).fullPivLu().kernel();
    res[i].normalize();
    if(res[i].dot(testVec) < 0){
      res[i] *= -1;
    }
  }
  return res;
};

void runHDGTransport(SimRun * thisRun,  HDGSolverType globType){
  //setup
  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order;
  thisRun->meshLocation += meshName + ".h5";
  std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/Transport/HDG/";
  if(globType == WEXPLICIT){
    writeDir += "WExp/";
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
  writeDir += meshName + "_dt-" + thisRun->timeStep;
  boost::filesystem::create_directory(writeDir);
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
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
  double invSqrt2 = 1.0/std::sqrt(2.0);
  for(int i = 0; i < vel.getLength(); i++){
    vel.getValues()->at(i) = invSqrt2;
  }
  Field tau(&myMesh, Face, nNodesPerFace, 2);
  Field anaSol(&myMesh, Node, 1, 1);
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
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGTransport transportMod(myMesh.getReferenceElement());
  double timeStep = std::stod(thisRun->timeStep);
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
  hdfio.setField("Solution", &sol);
  const std::set<int> * boundary = myMesh.getBoundaryFaces();
  std::vector<double> node(nodeDim);
  std::vector<int> cell(nNodesPerFace, 0);
  std::vector< std::vector<double> > nodes(nNodesPerFace, std::vector<double>(nodeDim, 0.0));
  double t = 0;
  std::vector<double> v(nodeDim, 0.0);
  std::vector<EVector> normals(nNodesPerFace, EVector::Zero(nodeDim));
  double projV = 0.0;
  for(int i = 0; i < myMesh.getNumberFaces(); i++){
    normals = calculateOutwardNormal(&myMesh, i);
    myMesh.getFace(i, &cell);
    for(int j = 0; j < nNodesPerFace; j++){
      vel.getValues(cell[j], &v);
      projV = normals[j].dot(EMap<EVector>(v.data(), v.size()));
      //tau.getValues()->at(2*(i*nNodesPerFace + j)) = 0.0;
      //tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = 0.0;
      //tau.getValues()->at(2*(i*nNodesPerFace + j)) = std::fabs(projV);
      //tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = std::fabs(projV);
      if(projV > 0){
        tau.getValues()->at(2*(i*nNodesPerFace + j)) = std::fabs(projV);
        tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = 0.0;
      } else {
        tau.getValues()->at(2*(i*nNodesPerFace + j)) = 0.0;
        tau.getValues()->at(2*(i*nNodesPerFace + j) + 1) = std::fabs(projV);
      }
    }
  }
  for(int i = 0; i < myMesh.getNumberCells(); i++){
    myMesh.getCell(i, &cell);
    for(int j = 0; j < nNodes; j++){
      myMesh.getPoint(cell[j], &node);
      vel.getValues(cell[j], &v);
      sol.getValues()->at(i*nNodes + j) = analyticalTransportHDG(t, node, v);
      for(int k = 0; k < nodeDim; k++){
        flux.getValues()->at((i*nNodes + j)*nodeDim + k) = analyticalTransportHDGGrad(t, node, v, k);
      }
    }
  }
  for(int i = 0; i < myMesh.getNumberFaces(); i++){
    myMesh.getFace(i, &cell);
    for(int j = 0; j < nNodesPerFace; j++){
      myMesh.getPoint(cell[j], &node);
      vel.getValues(cell[j], &v);
      trace.getValues()->at(i*nNodesPerFace + j) = analyticalTransportHDG(t, node, v);;
    }
  }
  hdfio.write(writeDir + "/res_0.h5");
  double timeEnd = 1.0;
  int nIters = timeEnd / timeStep;
  //int nIters = 2;
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
  int i = 0;
  ProgressBar pbar;
  pbar.setIterIndex(&i);
  pbar.setNumIterations(nIters);
  std::cout << "Simulation (d=" + thisRun->dim + ", h=" + thisRun->meshSize + ", p=" + thisRun->order + ", dt=" + thisRun->timeStep + ")" << std::endl;
  pbar.update();
  for(i = 0; i < nIters; i++){
    t += timeStep;
    //analytical sol
    for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
      myMesh.getPoint(iNode, &node);
      vel.getValues(iNode, &v);
      anaSol.getValues()->at(iNode) = analyticalTransportHDG(t, node, v);
    }
    for(auto it = boundary->begin(); it != boundary->end(); it++){
      myMesh.getFace(*it, &cell);
      myMesh.getSlicePoints(cell, &nodes);
      for(int j = 0; j < nNodesPerFace; j++){
        vel.getValues(cell[j], &v);
        dirichlet.getValues()->at((*it)*nNodesPerFace+j) = analyticalTransportHDG(t, nodes[j], v);
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
    for(int k = 0; k < ts.getNumStages(); k++){
      start = std::chrono::high_resolution_clock::now();
      mySolver.assemble();
      end = std::chrono::high_resolution_clock::now();
      thisRun->assembly += end - start;
      start = std::chrono::high_resolution_clock::now();
      mySolver.solve();
      ts.computeStage(&fieldMap);
      end = std::chrono::high_resolution_clock::now();
      thisRun->resolution += end - start;
    }
    start = std::chrono::high_resolution_clock::now();
    ts.computeSolution(&fieldMap);
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
      myMesh.getCell(k, &cell);
      for(int j = 0; j < nNodes; j++){
        residual.getValues()->at(k*nNodes + j) = anaSol.getValues()->at(cell[j]) - sol.getValues()->at(k*nNodes + j);
        sumRes += std::pow(residual.getValues()->at(k*nNodes + j), 2);
        sumAna += std::pow(anaSol.getValues()->at(cell[j]), 2);
      }
    }
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
      hdfio.write(writeDir + "/res_" + std::to_string(i+1) + ".h5");
    }
    end = std::chrono::high_resolution_clock::now();
    thisRun->post += end - start;
    pbar.update();
    if(thisRun->l2Err > 1.0){
      break;
    }
  }
};

TEST_CASE("Testing regression cases for Transport", "[regression][HDG][Transport]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  meshSizes["2"] = {"2e-1", "1e-1", "7e-2"};
  //meshSizes["2"] = {"1e-1"};
  std::vector<std::string> timeSteps = {"1e-2", "5e-3", "2e-3", "1e-3", "5e-4", "2e-4", "1e-4"};
  //std::vector<std::string> timeSteps = {"1e-3"};
  std::vector<std::string> orders = {"1", "2", "3"};
  //std::vector<std::string> orders = {"1"};
  std::vector<std::string> rkTypes = {"SSPRK3"};
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

  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/Transport/HDG/";
  //HDGSolverType globType = WEXPLICIT;
  HDGSolverType globType = IMPLICIT;
  std::string writeFile = "Breakdown.csv";
  if(globType == WEXPLICIT){
    writePath += "WExp/";
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
  std::ofstream f; f.open(writePath + writeFile);
  f << "dim,order,h,timeStep,linAlgErr,l2Err,dL2Err,runtime,setup,assembly,resolution,post\n";
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGTransport(&(*it), globType);
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    //CHECK(it->l2Err < 1e-1);
    f << it->dim << ",";
    f << it->order << ",";
    f << it->meshSize << ",";
    f << it->timeStep << ",";
    f << it->linAlgErr << ",";
    f << it->l2Err << ",";
    f << it->dL2Err << ",";
    f << it->runtime.count() << ",";
    f << it->setup.count() << ",";
    f << it->assembly.count() << ",";
    f << it->resolution.count() << ",";
    f << it->post.count() << "\n";
  }
  f.close();
};

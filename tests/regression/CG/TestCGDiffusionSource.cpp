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
#include "DiffusionSource.h"
#include "RungeKutta.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "CGSolver.h"
#include "TestUtils.h"

using namespace hfox;

double gaussianSrcCG(const std::vector<double> & v){
  double res = 0.0;
  double pre = 2.0/std::sqrt(M_PI);
  for(int i = 0; i < v.size(); i++){
    res -= pre*std::exp(-std::pow(v[i] - 0.5, 2));
  }
  return res;
};

double analyticalDiffSrcCG(const double t, const std::vector<double> & v){
  double res = 0.0;
  std::vector<std::complex<double> > w(v.size(), std::complex<double>(0.0, M_PI/2.0));
  double sqrtPi = std::sqrt(M_PI);
  std::complex<double> spaceArg(0.0, 0.0);
  double timePre = 0.0;
  double alpha = 0.0;
  for(int i = 0; i < v.size(); i++){
    spaceArg += w[i]*v[i];
    timePre += std::real(std::pow(w[i], 2));
    alpha = v[i] - 0.5;
    res += alpha*std::erf(alpha) + std::exp(-std::pow(alpha, 2))/sqrtPi;
  }
  res += std::exp(timePre * t)*std::real(std::exp(spaceArg));
  return res;
};

void runCGDiffSrc(SimRun * thisRun){
  //setup
  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order;
  thisRun->meshLocation += meshName + ".h5";
  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/DiffusionConvergence/CG/";
  std::string writeDir = writePath;
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
  int nNodes = myMesh.getReferenceElement()->getNumNodes();
  int nNodesPerFace = myMesh.getReferenceElement()->getFaceElement()->getNumNodes();
  Field dirichlet(&myMesh, Face, nNodesPerFace, 1);
  Field sol(&myMesh, Node, 1, 1);
  Field oldSol(&myMesh, Node, 1, 1);
  Field rk0(&myMesh, Node, 1, 1);
  Field rk1(&myMesh, Node, 1, 1);
  Field rk2(&myMesh, Node, 1, 1);
  Field rk3(&myMesh, Node, 1, 1);
  Field anaSol(&myMesh, Node, 1, 1);
  Field residual(&myMesh, Node, 1, 1);
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["OldSolution"] = &oldSol;
  fieldMap["RKStage_0"] = &rk0;
  fieldMap["RKStage_1"] = &rk1;
  fieldMap["RKStage_2"] = &rk2;
  fieldMap["RKStage_3"] = &rk3;
  fieldMap["Dirichlet"] = &dirichlet;
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  DiffusionSource diffMod(myMesh.getReferenceElement());
  RungeKutta ts(myMesh.getReferenceElement(), rkType);
  double timeStep = std::stod(thisRun->timeStep);
  ts.setTimeStep(timeStep);
  diffMod.setTimeScheme(&ts);
  PetscOpts myOpts;
  myOpts.maxits = 20000;
  myOpts.rtol = 1e-12;
  myOpts.verbose = false;
  PetscInterface petsciface(myOpts);
  CGSolver mySolver;
  mySolver.setVerbosity(0);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&diffMod);
  mySolver.setBoundaryModel(&dirMod);
  mySolver.initialize();
  mySolver.allocate();
  diffMod.setSourceFunction(gaussianSrcCG);
  hdfio.setField("Solution", &sol);
  hdfio.setField("OldSolution", &oldSol);
  hdfio.setField("RKStage_0", &rk0);
  hdfio.setField("RKStage_1", &rk1);
  hdfio.setField("RKStage_2", &rk2);
  const std::set<int> * boundary = myMesh.getBoundaryFaces();
  std::vector<double> node(myMesh.getNodeSpaceDimension());
  std::vector<int> cell(nNodesPerFace, 0);
  std::vector< std::vector<double> > nodes(nNodesPerFace, std::vector<double>(myMesh.getNodeSpaceDimension(), 0.0));
  double t = 0;
  for(int i = 0; i < myMesh.getNumberPoints(); i++){
    myMesh.getPoint(i, &node);
    sol.getValues()->at(i) = analyticalDiffSrcCG(t, node);
  }
  hdfio.write(writeDir + "/res_0.h5");
  double timeEnd = 1.0;
  int nIters = timeEnd / timeStep;
  //int nIters = 5;
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
      anaSol.getValues()->at(iNode) = analyticalDiffSrcCG(t, node);
    }
    for(auto it = boundary->begin(); it != boundary->end(); it++){
      myMesh.getFace(*it, &cell);
      myMesh.getSlicePoints(cell, &nodes);
      for(int j = 0; j < nNodesPerFace; j++){
        dirichlet.getValues()->at((*it)*nNodesPerFace+j) = analyticalDiffSrcCG(t, nodes[j]);
      }
    }
    //copy sol into oldsol
    for(int iN = 0; iN < myMesh.getNumberPoints(); iN++){
      oldSol.getValues()->at(iN) = sol.getValues()->at(iN);
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
    for(int k = 0; k < myMesh.getNumberPoints(); k++){
      residual.getValues()->at(k) = anaSol.getValues()->at(k) - sol.getValues()->at(k);
      sumRes += std::pow(residual.getValues()->at(k), 2);
      sumAna += std::pow(anaSol.getValues()->at(k), 2); 
    }
    double l2Err = std::sqrt(TestUtils::l2ProjectionNodeField(&residual, &residual, &myMesh));
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
      end = std::chrono::high_resolution_clock::now();
      thisRun->post += end - start;
    }
    pbar.update();
  }
};

TEST_CASE("Testing regression cases for DiffusionSource", "[regression][CG][DiffusionSource]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  meshSizes["2"] = {"2e-1", "1e-1", "7e-2"};
  //meshSizes["2"] = {"1e-1"};
  std::vector<std::string> timeSteps = {"2e-1", "1e-1", "5e-2", "1e-2", "5e-3", "2e-3", "1e-3", "5e-4", "2e-4"};
  //std::vector<std::string> timeSteps = {"2e-1", "1e-1", "5e-2", "1e-2", "5e-3", "2e-3"};
  std::vector<std::string> orders = {"1", "2", "3"};
  //std::vector<std::string> orders = {"3"};
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

  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/DiffusionConvergence/CG/";
  std::string writeFile = "Breakdown.csv";
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
    runCGDiffSrc(&(*it));
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    //CHECK(it->l2Err < 1e-2);
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

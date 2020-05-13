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
#include "HDGDiffusionSource.h"
#include "Euler.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "TestUtils.h"

using namespace hfox;

double gaussianSrc(const std::vector<double> & v){
  double res = 0.0;
  double pre = 2.0/std::sqrt(M_PI);
  for(int i = 0; i < v.size(); i++){
    res -= pre*std::exp(-std::pow(v[i] - 0.5, 2));
  }
  return res;
};

double analyticalDiffSrc(const double t, const std::vector<double> & v){
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

void runHDGDiffSrc(SimRun * thisRun, bool isExplicit, HDGSolverType globType){
  //setup
  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  thisRun->meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/";
  std::string meshName = "regression_dim-" + thisRun->dim + "_h-" + thisRun->meshSize;
  meshName += "_ord-" + thisRun->order;
  thisRun->meshLocation += meshName + ".h5";
  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/DiffusionConvergence/";
  std::string writeDir = writePath;
  if(!isExplicit){
    switch(globType){
      case IMPLICIT:{writeDir += "ImpImp/";break;}
      case WEXPLICIT:{writeDir += "WExpImp/";break;}
    }
  } else{
    switch(globType){
      case IMPLICIT:{writeDir += "ImpExp/";break;}
      case WEXPLICIT:{writeDir += "WExpExp/";break;}
    }
  }
  writeDir += meshName + "_dt-" + thisRun->timeStep;
  boost::filesystem::create_directory(writeDir);
  Mesh myMesh(std::stoi(thisRun->dim), std::stoi(thisRun->order), "simplex");
  HDF5Io hdfio(&myMesh);
  hdfio.load(thisRun->meshLocation);
  int nNodes = myMesh.getReferenceElement()->getNumNodes();
  int nNodesPerFace = myMesh.getReferenceElement()->getFaceElement()->getNumNodes();
  Field dirichlet(&myMesh, Face, nNodesPerFace, 1);
  Field sol(&myMesh, Cell, nNodes, 1);
  Field flux(&myMesh, Cell, nNodes, myMesh.getNodeSpaceDimension());
  Field trace(&myMesh, Face, nNodesPerFace, 1);
  Field tau(&myMesh, Face, nNodesPerFace, 1);
  std::fill(tau.getValues()->begin(), tau.getValues()->end(), 1.0);
  Field anaSol(&myMesh, Node, 1, 1);
  Field residual(&myMesh, Cell, nNodes, 1);
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["Flux"] = &flux;
  fieldMap["Trace"] = &trace;
  fieldMap["Tau"] = &tau;
  fieldMap["Dirichlet"] = &dirichlet;
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  HDGDiffusionSource diffMod(myMesh.getReferenceElement());
  Euler ts(myMesh.getReferenceElement(), isExplicit);
  double timeStep = std::stod(thisRun->timeStep);
  ts.setTimeStep(timeStep);
  diffMod.setTimeScheme(&ts);
  PetscOpts myOpts;
  myOpts.maxits = 20000;
  myOpts.rtol = 1e-6;
  myOpts.verbose = true;
  PetscInterface petsciface(myOpts);
  HDGSolverOpts solveOpts;
  solveOpts.type = globType;
  solveOpts.verbosity = true;
  HDGSolver mySolver;
  mySolver.setOptions(solveOpts);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&diffMod);
  mySolver.setBoundaryModel(&dirMod);
  mySolver.initialize();
  mySolver.allocate();
  diffMod.setSourceFunction(gaussianSrc);
  hdfio.setField("Solution", &sol);
  const std::set<int> * boundary = myMesh.getBoundaryFaces();
  std::vector<double> node(myMesh.getNodeSpaceDimension());
  std::vector<int> cell(nNodesPerFace, 0);
  std::vector< std::vector<double> > nodes(nNodesPerFace, std::vector<double>(myMesh.getNodeSpaceDimension(), 0.0));
  double t = 0;
  for(int i = 0; i < myMesh.getNumberCells(); i++){
    myMesh.getCell(i, &cell);
    for(int j = 0; j < nNodes; j++){
      myMesh.getPoint(cell[j], &node);
      sol.getValues()->at(i*nNodes + j) = analyticalDiffSrc(t, node);
    }
  }
  hdfio.write(writeDir + "/res_0.h5");
  double timeEnd = 5*timeStep;
  //double timeEnd = 1.0;
  int nIters = timeEnd / timeStep;
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
  for(int i = 0; i < nIters; i++){
    t += timeStep;
    //analytical sol
    for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
      myMesh.getPoint(iNode, &node);
      anaSol.getValues()->at(iNode) = analyticalDiffSrc(t, node);
    }
    for(auto it = boundary->begin(); it != boundary->end(); it++){
      myMesh.getFace(*it, &cell);
      myMesh.getSlicePoints(cell, &nodes);
      for(int j = 0; j < nNodesPerFace; j++){
        dirichlet.getValues()->at((*it)*nNodesPerFace+j) = analyticalDiffSrc(t, nodes[j]);
      }
    }
    start = std::chrono::high_resolution_clock::now();
    mySolver.assemble();
    end = std::chrono::high_resolution_clock::now();
    thisRun->assembly += end - start;
    start = std::chrono::high_resolution_clock::now();
    mySolver.solve();
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
    std::cout << "l2Err at time " << t << ": " << l2Err << std::endl;
    if(i != (nIters-1)){
      thisRun->linAlgErr += linAlgErr*timeStep;
      thisRun->l2Err += l2Err*timeStep;
      thisRun->dL2Err += dL2Err*timeStep;
    } else{
      thisRun->linAlgErr += linAlgErr*timeStep/2;
      thisRun->l2Err += l2Err*timeStep/2;
      thisRun->dL2Err += dL2Err*timeStep/2;
    }
    //double quot = t/(5e-3);
    double quot = 0.0;
    double rem = quot - ((int)quot);
    std::cout << "rem: " << rem << std::endl;
    if(rem < timeStep/(5e-3)){
      hdfio.write(writeDir + "/res_" + std::to_string(i+1) + ".h5");
      end = std::chrono::high_resolution_clock::now();
      thisRun->post += end - start;
    }
  }
};

TEST_CASE("Testing regression cases for HDGDiffusionSource", "[regression][HDG][DiffusionSource]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  //meshSizes["2"] = {"2e-1", "1e-1", "7e-2"};
  meshSizes["2"] = {"5e-2"};
  //std::vector<std::string> timeSteps = {"2e-1", "1e-1", "5e-2", "1e-2", "5e-3"};
  std::vector<std::string> timeSteps = {"1e-11"};
  //std::vector<std::string> orders = {"1", "2", "3", "4", "5"};
  std::vector<std::string> orders = {"3"};
  std::vector<SimRun> simRuns;
  for(auto it = meshSizes.begin(); it != meshSizes.end(); it++){
    for(auto itMs = it->second.begin(); itMs != it->second.end(); itMs++){
      for(int o = 0; o < orders.size(); o++){
        for(int dt = 0; dt < timeSteps.size(); dt++){
          SimRun thisRun;
          thisRun.dim = it->first;
          thisRun.meshSize = *itMs;
          thisRun.order = orders[o];
          thisRun.timeStep = timeSteps[dt];
          simRuns.push_back(thisRun);
        }
      }
    }
  }

  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/DiffusionConvergence/";
  bool isExplicit = 1;
  HDGSolverType globType = WEXPLICIT;
  //HDGSolverType globType = IMPLICIT;
  std::string writeFile = "Breakdown.csv";
  if(!isExplicit){
    switch(globType){
      case IMPLICIT:{writePath += "ImpImp/";break;}
      case WEXPLICIT:{writePath += "WExpImp/";break;}
    }
  } else{
    switch(globType){
      case IMPLICIT:{writePath += "ImpExp/";break;}
      case WEXPLICIT:{writePath += "WExpExp/";break;}
    }
  }
  std::ofstream f; f.open(writePath + writeFile);
  f << "dim,order,h,timeStep,linAlgErr,l2Err,dL2Err,runtime,setup,assembly,resolution,post\n";
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runHDGDiffSrc(&(*it), isExplicit, globType);
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

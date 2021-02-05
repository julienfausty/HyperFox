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
#include "HDGBohmModel.h"
#include "IntegratedDirichletModel.h"
#include "PetscInterface.h"
#include "HDGSolver.h"
#include "NonLinearWrapper.h"
#include "TestUtils.h"

using namespace hfox;

using namespace nGamma;

double initDensity(std::vector<double> & node){
  return 1.0;
};//initDensity

double initGamma(std::vector<double> & node){
  //return (-1.0 + node[0]*2.0);
  return 0.0;
  //return 2.0*(node[0]-0.5) - std::sin(M_PI*node[0]);
};//initGamma

double numHeavyside(double x){
  double delta = 0.001;
  return 0.5*(1.0 + std::tanh(x/delta));
};

double numDerivHeavyside(double x){
  double delta = 0.001;
  return (0.5/delta)*(1.0 - std::pow(std::tanh(x/delta), 2));
};

double tanhTransferF(double x, double y){
  double Hx = numHeavyside(x);
  double Hy = numHeavyside(y);
  double Hmx = numHeavyside(-x);
  double Hmy = numHeavyside(-y);
  double res = 0.0;
  double tol = 1e-5;
  res += (Hx*Hmy + Hmx*Hy);
  res += Hx*Hy*numHeavyside(y-x);
  res += Hmx*Hmy*numHeavyside(x-y);
  res *= (1.0 - numHeavyside(std::abs(tol) - std::abs(y)));
  return res;
};

std::vector<double> tanhDerivTransferF(double x, double y){
  std::vector<double> res(2, 0);
  double Hx = numHeavyside(x);
  double Hy = numHeavyside(y);
  double Hmx = numHeavyside(-x);
  double Hmy = numHeavyside(-y);
  double dHx = numDerivHeavyside(x);
  double dHy = numDerivHeavyside(y);
  double dHmx = numDerivHeavyside(-x);
  double dHmy = numDerivHeavyside(-y);
  double tol = 1e-5;
  double Htol = numHeavyside(std::abs(tol) - std::abs(y));
  double dHtol = numDerivHeavyside(std::abs(tol) - std::abs(y));
  double transfer = (Hx*Hmy + Hmx*Hy + Hx*Hy*numHeavyside(y-x) +  Hmx*Hmy*numHeavyside(x-y));
  res[0] += dHx*Hmy - dHmx*Hy;
  res[0] += dHx*Hy*numHeavyside(y-x) - Hx*Hy*numDerivHeavyside(y-x);
  res[0] += -dHmx*Hmy*numHeavyside(x-y) + Hmx*Hmy*numDerivHeavyside(x-y); 
  res[0] *= (1-Htol);
  res[1] += -Hx*dHmy + Hmx*dHy;
  res[1] += Hx*dHy*numHeavyside(y-x) + Hx*Hy*numDerivHeavyside(y-x);
  res[1] += -Hmx*dHmy*numHeavyside(x-y) - Hmx*Hmy*numDerivHeavyside(x-y);
  res[1] *= (1-Htol);
  res[1] += transfer*dHtol;
  return res;
};

double heavysideTransfer(double x, double y){
  double res = 0.0;
  double tol = 1e-5;
  bool xsup = x > 0;
  bool ysup = y > 0;
  bool ytol = (std::abs(y) > tol);
  if(ytol){
    if((xsup != ysup)){
      res = 1.0;
    } else {
      if(std::abs(x) - std::abs(y) < 0.0){
        res = 1.0;
      }
    }
  }
  return res;
};

std::vector<double> heavysideDerivTransfer(double x, double y){
  std::vector<double> res(2, 0);
  return res;
};

double neumannTransfer(double x, double y){
  double res = 0.0;
  return res;
};

std::vector<double> neumannDerivTransfer(double x, double y){
  std::vector<double> res(2, 0);
  return res;
};

double dirichletTransfer(double x, double y){
  double res = 1.0;
  return res;
};

std::vector<double> dirichletDerivTransfer(double x, double y){
  std::vector<double> res(2, 0);
  return res;
};

double constantSrc(std::vector<double> x, int i){
  double res = 0.0;
  if(i == 0){
    res = 2.0;
  }
  return res;
};//gaussianSrcnGamma

double gaussianSrcnGamma(std::vector<double> x, int i){
  double res = 0.0;
  if(i == 0){
    double sigma = 0.01;
    res = std::exp(-std::pow(x[0]-0.5, 2)/sigma);
  }
  return res;
};//gaussianSrcnGamma

void calculateBoundaryNormals(Mesh * myMesh, int iFace, std::vector<double> * normals){
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesPEl = refEl->getNumNodes();
  int nNodesPFc = fEl->getNumNodes();
  int dimNode = myMesh->getNodeSpaceDimension();
  std::vector<int> cell(nNodesPEl, 0);
  std::vector<int> face(nNodesPFc, 0);
  std::vector<int> face2Cell(1, 0);
  std::vector< std::vector<double> > faceNodes(nNodesPFc, std::vector<double>(dimNode, 0.0));
  EMatrix faceNodeMat(dimNode, nNodesPFc);
  EMatrix shapeDerivMat(nNodesPFc, dimNode-1);
  std::vector<EMatrix> jacobians(nNodesPFc, EMatrix::Zero(dimNode, dimNode-1));
  normals->resize(dimNode*nNodesPFc, 0.0);
  myMesh->getFace(iFace, &face);
  myMesh->getFace2Cell(iFace, &face2Cell);
  myMesh->getCell(face2Cell[0], &cell);
  myMesh->getSlicePoints(face, &faceNodes);
  for(int iN = 0; iN < nNodesPFc; iN++){
    faceNodeMat.col(iN) = EMap<EVector>(faceNodes.at(iN).data(), dimNode);
  }
  for(int iN = 0; iN < nNodesPFc; iN++){
    for(int jN = 0; jN < nNodesPFc; jN++){
      shapeDerivMat.row(jN) = EMap<const EVector>(fEl->getDerivShapeFunctions()->at(iN)[jN].data(), dimNode-1);
    }
    jacobians[iN] = faceNodeMat*shapeDerivMat;
  }
  EVector vn(dimNode);
  for(int iN = 0; iN < nNodesPEl; iN++){
    std::vector<int>::iterator itv = std::find(face.begin(), face.end(), cell[iN]);
    if(itv == face.end()){
      int locInd = myMesh->getPartitioner()->global2LocalNode(cell[iN]);
      if(locInd != -1){
        std::vector<double> node(dimNode, 0.0);
        myMesh->getPoint(locInd, &node);
        vn = EMap<EVector>(node.data(), node.size()) - EMap<EVector>(faceNodes[0].data(), faceNodes[0].size());
        break;
      }
    }
  }
  for(int iN = 0; iN < nNodesPFc; iN++){
    EMatrix temp = jacobians[iN].transpose().fullPivLu().kernel();
    EMap<EVector> normal(normals->data() + iN*dimNode, dimNode);
    normal = temp.col(0);
    normal.normalize();
    double prod = vn.dot(normal);
    if(prod > 0){
      normal *= -1;
    }
  }
};

void runnGammaBohmSimulation(SimRun * thisRun, HDGSolverType globType){
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
  Field sol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field buffSol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field oldSol(&myMesh, Cell, nNodesPerEl, nodeDim);
  Field flux(&myMesh, Cell, nNodesPerEl, std::pow(nodeDim, 2));
  Field oldFlux(&myMesh, Cell, nNodesPerEl, std::pow(nodeDim, 2));
  Field trace(&myMesh, Face, nNodesPerFace, nodeDim);
  Field oldTrace(&myMesh, Face, nNodesPerFace, nodeDim);
  Field partition(&myMesh, Node, 1, 1);
  Field extN(&myMesh, Face, nNodesPerFace, nodeDim); 
  Field dirichlet(&myMesh, Face, nNodesPerFace, nodeDim); 
  Field b(&myMesh, Node, 1, nodeDim);
  Field D(&myMesh, Node, 1, 1);
  Field G(&myMesh, Node, 1, 1);
  Field c(&myMesh, Node, 1, 1);
  Field tau(&myMesh, Face, nNodesPerFace, std::pow(nodeDim, 2)*2);
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
  fieldMap["b"] = &b;
  fieldMap["D"] = &D;
  fieldMap["G"] = &G;
  fieldMap["SoundVelocity"] = &c;
  fieldMap["ExteriorNormals"] = &extN;
  fieldMap["Dirichlet"] = &dirichlet;
  //initialize models, solvers, etc.
  IntegratedDirichletModel dMod(myMesh.getReferenceElement()->getFaceElement());
  std::set<int> dfaces;
  HDGBohmModel bMod(myMesh.getReferenceElement()->getFaceElement());
  std::set<int> bfaces;
  //bMod.setTransferFunction(tanhTransferF, tanhDerivTransferF);
  bMod.setTransferFunction(heavysideTransfer, heavysideDerivTransfer);
  //bMod.setTransferFunction(neumannTransfer, neumannDerivTransfer);
  //bMod.setTransferFunction(dirichletTransfer, dirichletDerivTransfer);
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
  mySolver.setBoundaryCondition(&dMod, &dfaces);
  mySolver.setBoundaryCondition(&bMod, &bfaces);
  NonLinearWrapper wrapper;
  wrapper.setVerbosity(false);
  int maxNRIter = 15;
  wrapper.setMaxIterations(maxNRIter);
  wrapper.setResidualTolerance(1e-6);
  wrapper.setSolutionFields(&sol, &buffSol);
  wrapper.setSolver(&mySolver);
  //setup outputs
  std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/BohmResults/";
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
  hdfio.setField("Trace", &trace);
  hdfio.setField("Partition", &partition);
  //define scalars
  double diffCoeff = 0.01;//diffusive coeff
  double timeStep = std::stod(thisRun->timeStep);
  double carLen = 1.0;//Imp
  double t = 0;
  double timeEnd = 1.0;
  double cs = 1.0;
  int nIters = timeEnd/timeStep;
  nGammaParams locParams;
  locParams.soundSpeed = cs;
  transportMod.setParams(locParams);
  //create buffers
  std::vector<int> cell;
  std::vector<double> node;
  std::vector<int> face2Cell;
  std::vector<int> ibuffer;
  std::vector<double> dbuffer;
  int ibuff;
  //initialize field values
  for(int i = 0; i < myMesh.getNumberCells(); i++){
    myMesh.getCell(i, &cell);
    for(int j = 0; j < cell.size(); j++){
      myMesh.getPoint(cell[j], &node);
      sol.getValues()->at((i*nNodesPerEl + j)*nodeDim) = initDensity(node);//density
      sol.getValues()->at((i*nNodesPerEl + j)*nodeDim + 1) = initGamma(node);//parallel momentum
    }
  }
  std::copy(sol.getValues()->begin(), sol.getValues()->end(), buffSol.getValues()->begin());
  std::fill(D.getValues()->begin(), D.getValues()->end(), diffCoeff);
  std::fill(G.getValues()->begin(), G.getValues()->end(), diffCoeff);
  std::fill(c.getValues()->begin(), c.getValues()->end(), cs);
  for(int iNode = 0; iNode < myMesh.getNumberPoints(); iNode++){
    b.getValues()->at(iNode*nodeDim) = 0.1;
    b.getValues()->at(iNode*nodeDim + 1) = 0.0;
  }
  for(int iFace = 0; iFace < myMesh.getNumberFaces(); iFace++){
    myMesh.getFace(iFace, &cell);
    for(int iN = 0; iN < nNodesPerFace; iN++){
      for(int iCell = 0; iCell < 2; iCell++){
        //EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          //= EMatrix::Identity(nodeDim, nodeDim)*(diffCoeff/carLen)*2.0;
        EMap<EMatrix>(tau.getValues()->data() + ((iFace*nNodesPerFace + iN)*2 + iCell)*nodeDim*nodeDim, nodeDim, nodeDim)
          = EMatrix::Identity(nodeDim, nodeDim);
      }
    }
  }
  for(std::set<int>::const_iterator itFc = myMesh.getBoundaryFaces()->begin(); itFc != myMesh.getBoundaryFaces()->end(); itFc++){
    myMesh.getFace(*itFc, &ibuffer);
    calculateBoundaryNormals(&myMesh, *itFc, &dbuffer);
    std::copy(dbuffer.begin(), dbuffer.end(), extN.getValues()->begin() + (*itFc)*nNodesPerFace*nodeDim);
    for(int k = 0; k < nNodesPerFace; k++){
      myMesh.getPoint(ibuffer[k], &node);
      dirichlet.getValues()->at(((*itFc)*nNodesPerFace + k)*nodeDim) = initDensity(node);
      dirichlet.getValues()->at(((*itFc)*nNodesPerFace + k)*nodeDim + 1) = initGamma(node);
      trace.getValues()->at(((*itFc)*nNodesPerFace + k)*nodeDim) = initDensity(node);
      trace.getValues()->at(((*itFc)*nNodesPerFace + k)*nodeDim + 1) = initGamma(node);
    }
  }
  //create field list
  std::vector<Field*> fieldList;
  for(auto itMap = fieldMap.begin(); itMap != fieldMap.end(); itMap++){
    fieldList.push_back(itMap->second);
  }
  fieldList.push_back(&partition);
  //partition the mesh and fields
  zPart.setFields(fieldList);
  zPart.computePartition();
  zPart.update();
  //fill partition field
  for(int iN = 0; iN < myMesh.getNumberPoints(); iN++){
    partition.getValues()->at(iN) = zPart.getRank();
  }
  //setup boundary faces
  for(std::set<int>::const_iterator itFc = myMesh.getBoundaryFaces()->begin(); itFc != myMesh.getBoundaryFaces()->end(); itFc++){
    int ibuff = zPart.global2LocalFace(*itFc);
    myMesh.getFace(ibuff, &ibuffer);
    bool isDFace = 1;
    for(int k = 0; k < nNodesPerFace; k++){
      ibuff = zPart.global2LocalNode(ibuffer[k]);
      if(ibuff != -1){
        myMesh.getPoint(ibuff, &node);
      } else {
        myMesh.getGhostPoint(ibuffer[k], &node);
      }
      isDFace = (node[0] != 1.0);
      if(isDFace){
        break;
      }
    }
    isDFace = 0; //full Bohm
    //isDFace = 1; //full Dirichlet
    if(isDFace){
      dfaces.insert(*itFc);
    } else {
      bfaces.insert(*itFc);
    }
  }
  //first output
  //hdfio.write(writeDir + "/res_0.h5");
  //allocating and last set ups
  ts.setTimeStep(timeStep);
  transportMod.setTimeScheme(&ts);
  mySolver.initialize();
  mySolver.allocate();
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
    //transportMod.setSourceFunction(gaussianSrcnGamma);
    transportMod.setSourceFunction(constantSrc);
    //compute necessary fields
    std::copy(sol.getValues()->begin(), sol.getValues()->end(), buffSol.getValues()->begin());
    std::copy(sol.getValues()->begin(), sol.getValues()->end(), oldSol.getValues()->begin());
    std::copy(flux.getValues()->begin(), flux.getValues()->end(), oldFlux.getValues()->begin());
    std::copy(trace.getValues()->begin(), trace.getValues()->end(), oldTrace.getValues()->begin());
    zPart.updateSharedInformation();
    //solve non-linear problem
    start = std::chrono::high_resolution_clock::now();
    for(int k = 0; k < ts.getNumStages(); k++){
      wrapper.solve();//Imp
      //mySolver.assemble();//Exp
      //mySolver.solve();//Exp
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
    if(i != (nIters-1)){
      thisRun->linAlgErr += linAlgErr*timeStep;
    } else{
      thisRun->linAlgErr += linAlgErr*timeStep/2;
    }
    //double quot = t/(1e-3);
    //double rem = quot - ((int)quot);
    //if(rem < timeStep/(1e-3)){
      //hdfio.write(writeDir + "/res_" + std::to_string(i+1) + ".h5");
    //}
    end = std::chrono::high_resolution_clock::now();
    thisRun->post += end - start;
    //pbar.update();
  }
};


TEST_CASE("Testing regression cases for the HDGBohmModel", "[regression][HDG][BohmModel]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  meshSizes["2"] = {"1e-1"};
  //meshSizes["2"] = {"7e-2"};
  //std::vector<std::string> timeSteps = {"1e-2", "5e-3", "2e-3", "1e-3", "5e-4", "2e-4", "1e-4", "5e-5", "2e-5"};
  std::vector<std::string> timeSteps = {"1e-2"};
  //std::vector<std::string> timeSteps = {"1e-3"};
  std::vector<std::string> orders = {"3"};
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

  std::string writePath =  "/home/julien/workspace/M2P2/Postprocess/results/BohmResults/";
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
    runnGammaBohmSimulation(&(*it), globType);
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    it->runtime = end - start;
    CHECK(it->linAlgErr < 1);
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

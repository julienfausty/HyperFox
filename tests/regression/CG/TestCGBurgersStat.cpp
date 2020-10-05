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
#include "BurgersModel.h"
#include "DirichletModel.h"
#include "PetscInterface.h"
#include "CGSolver.h"
#include "NonLinearWrapper.h"
#include "TestUtils.h"

using namespace hfox;

namespace cg{

void analyticalBurgersStat(const std::vector<double> & x, std::vector<EMatrix> cMat, std::vector<double> * sol){
  sol->resize(x.size(), 0.0);
  double buff;
  for(int i = 0; i < cMat.size(); i++){
    sol->at(i) = 0.0;
    for(int j = 0; j < cMat[i].cols(); j++){
      buff = 1.0;
      for(int d = 0; d < x.size(); d++){
        buff *= std::pow(x[d], cMat[i](d, j));
      }
      buff *= cMat[i](x.size(), j);
      sol->at(i) += buff;
    }
  }
};

double sourceBurgersStat(const std::vector<double> & x, int i, std::vector<EMatrix> cMat, double D){
  double res = 0.0;
  double cBuff, dBuff;
  double ccBuff;
  for(int j = 0; j < x.size(); j++){
    for(int I = 0; I < cMat[i].cols(); I++){
      ccBuff = 0.0;
      for(int J = 0; J < cMat[j].cols(); J++){
        cBuff = 1.0;
        for(int d = 0; d < x.size(); d++){
          cBuff *= std::pow(x[d], cMat[j](d, J));
          if((d == j)){
            if(cMat[i](j, I) > 0){
              cBuff *= std::pow(x[d], cMat[i](j, I) - 1);
            } else {
              cBuff *= 0.0;
              break;
            }
          } else {
            cBuff *= std::pow(x[d], cMat[i](d, I));
          }
        }
        cBuff *= cMat[j](x.size(), J);
        ccBuff += cBuff;
      }
      dBuff = 1.0;
      for(int d = 0; d < x.size(); d++){
        if((d == j)){
          if(cMat[i](j, I) > 1){
            dBuff *= D * (cMat[i](j, I) - 1) * std::pow(x[d], cMat[i](j, I) - 2);
          } else {
            dBuff *= 0.0;
            break;
          }
        } else {
          dBuff *= std::pow(x[d], cMat[i](d, I));
        }
      }
      res += cMat[i](j, I)*cMat[i](x.size(), I)*(ccBuff - dBuff);
    }
  }
  return res;
};

};//cg

void runCGBurgersStat(SimRun * thisRun){
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
  Field sol(&myMesh, Node, 1, nodeDim);
  Field buffSol(&myMesh, Node, 1, nodeDim);
  Field partition(&myMesh, Node, 1, 1);
  Field diffCoeff(&myMesh, Node, 1, 1);
  Field anaSol(&myMesh, Node, 1, nodeDim);
  Field residual(&myMesh, Node, 1, 1);
  //create fieldMap
  std::map<std::string, Field*> fieldMap;
  fieldMap["Solution"] = &sol;
  fieldMap["BufferSolution"] = &buffSol;
  fieldMap["Dirichlet"] = &dirichlet;
  fieldMap["DiffusionTensor"] = &diffCoeff;
  //initialize models, solvers, etc.
  DirichletModel dirMod(myMesh.getReferenceElement()->getFaceElement());
  BurgersModel transportMod(myMesh.getReferenceElement());
  PetscOpts myOpts;
  myOpts.maxits = 6000;
  myOpts.rtol = 1e-16;
  myOpts.verbose = false;
  PetscInterface petsciface(myOpts);
  CGSolver mySolver;
  mySolver.setVerbosity(false);
  mySolver.setMesh(&myMesh);
  mySolver.setFieldMap(&fieldMap);
  mySolver.setLinSystem(&petsciface);
  mySolver.setModel(&transportMod);
  mySolver.setBoundaryModel(&dirMod);
  NonLinearWrapper wrapper;
  int maxIters = 50;
  wrapper.setResidualTolerance(1e-3);
  wrapper.setMaxIterations(maxIters);
  wrapper.setSolutionFields(&sol, &buffSol);
  wrapper.setSolver(&mySolver);
  //setup outputs
  //std::string writeDir = "/home/jfausty/workspace/Postprocess/results/BurgersStat/CG/";
  std::string writeDir = "/home/julien/workspace/M2P2/Postprocess/results/BurgersStat/CG/";
  writeDir += meshName;
  if(zPart.getRank() == 0){
    boost::filesystem::create_directory(writeDir);
  }
  hdfio.setField("Solution", &sol);
  hdfio.setField("Analytical", &anaSol);
  hdfio.setField("Residual", &residual);
  hdfio.setField("Partition", &partition);
  //define scalars
  double D = 1.0;//diffusive coeff 
  int dimGrad = std::pow(nodeDim, 2);
  std::vector<EMatrix> cMat(nodeDim, EMatrix::Zero(nodeDim+1, 1));
  cMat[0](2, 0) = 1.0;
  cMat[0](0, 0) = 5.0;
  //create buffers
  std::vector<int> cell;
  std::vector<double> node;
  std::vector<double> dbuffer;
  //initialize field values
  for(int i = 0; i < myMesh.getNumberPoints(); i++){
    myMesh.getPoint(i, &node);
    cg::analyticalBurgersStat(node, cMat, &dbuffer);
    for(int k = 0; k < dbuffer.size(); k++){
      anaSol.getValues()->at(i*nodeDim + k) = dbuffer[k];
    }
  }
  std::fill(diffCoeff.getValues()->begin(), diffCoeff.getValues()->end(), D);
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
  transportMod.setSourceFunction([cMat, D](const std::vector<double> & x, int i){return cg::sourceBurgersStat(x, i, cMat, D);});
  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
  thisRun->setup = end - start;
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
        cg::analyticalBurgersStat(node, cMat, &dbuffer);
        for(int dof = 0; dof < nodeDim; dof++){
          dirichlet.getValues()->at((locFace * nNodesPerFace + k)*nodeDim + dof) = dbuffer[dof];
        }
      }
    }
  }
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
  for(int k = 0; k < myMesh.getNumberPoints(); k++){
    residual.getValues()->at(k) = 0.0;
    for(int dof = 0; dof < nodeDim; dof++){
      residual.getValues()->at(k) += std::pow(anaSol.getValues()->at(k*nodeDim + dof) - sol.getValues()->at(k*nodeDim + dof), 2);
    }
    residual.getValues()->at(k) = std::sqrt(residual.getValues()->at(k));
    sumRes += std::pow(residual.getValues()->at(k), 2);
    sumAna += std::pow(anaSol.getValues()->at(k), 2); 
  }
  double l2Err = std::sqrt(TestUtils::l2ProjectionNodeField(&residual, &residual, &myMesh));
  double dL2Err = std::sqrt(sumRes/sumAna);
  //std::cout << "l2Err at time " << t << ": " << l2Err << std::endl;
  thisRun->linAlgErr = linAlgErr;
  thisRun->l2Err = l2Err;
  thisRun->dL2Err = dL2Err;
  hdfio.write(writeDir + "/res_" + std::to_string(maxIters) + ".h5");
  end = std::chrono::high_resolution_clock::now();
  thisRun->post += end - start;

};



TEST_CASE("Testing stationary regression cases for the BurgersModel", "[regression][CG][BurgersModelStat]"){
  std::map<std::string, std::vector<std::string> > meshSizes;
  //meshSizes["3"] = {"3e-1", "2e-1", "1e-1"};
  //meshSizes["2"] = {"3e-1", "2e-1", "1e-1", "7e-2", "5e-2"};
  //meshSizes["3"] = {"3e-1"};
  //meshSizes["2"] = {"7e-2"};
  meshSizes["2"] = {"2e-1", "1e-1", "7e-2", "5e-2", "2e-2"};
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

  //std::string writePath = "/home/jfausty/workspace/Postprocess/results/Burgers/CG/";
  std::string writePath = "/home/julien/workspace/M2P2/Postprocess/results/BurgersStat/CG/";
  std::string writeFile = "Breakdown.csv"; 
  int nParts;
  MPI_Comm_size(MPI_COMM_WORLD, &nParts);
  //writePath += std::to_string(nParts) + "/";
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ofstream f;
  if(rank == 0){
    f.open(writePath + writeFile);
    f << "dim,order,h,linAlgErr,l2Err,dL2Err,avgRuntime,maxRuntime,minRuntime\n" << std::flush;
  }
  for(auto it = simRuns.begin(); it != simRuns.end(); it++){
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    runCGBurgersStat(&(*it));
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

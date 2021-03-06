#include "NonLinearWrapper.h"

namespace hfox{

NonLinearWrapper::NonLinearWrapper() : mySolver(NULL), currentSolution(NULL), previousSolution(NULL), resTol(1e-6), maxIters(1000), dampening(0.0){
  setResidualComputer(NonLinearWrapper::vanillaResidualComputer);
  setLinearizedSolver(NonLinearWrapper::vanillaLinearizedSolver);
};//constructor

NonLinearWrapper::~NonLinearWrapper(){
};//destructor

double NonLinearWrapper::vanillaResidualComputer(Field * currentSol, Field * prevSol){
  double res = 0;
  double diffSum= 0;
  double refSum = 0;
  const std::vector<double> * curPtr, * prevPtr;
  curPtr = currentSol->getValues();
  prevPtr = prevSol->getValues();
  for(int i = 0; i < curPtr->size(); i++){
    diffSum += std::pow(curPtr->at(i) - prevPtr->at(i), 2.0);
    refSum += std::pow(prevPtr->at(i), 2.0);
  }
  double globDiffSum = 0;
  double globRefSum = 0;
  MPI_Allreduce(&diffSum, &globDiffSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&refSum, &globRefSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(globRefSum != 0){
    res = std::sqrt(globDiffSum/globRefSum);
  } else {
    res = std::sqrt(globDiffSum);
  }
  return res;
};//vanillaResidualComputer

void NonLinearWrapper::vanillaLinearizedSolver(Solver * solver){
  solver->assemble();
  solver->solve();
};//vanillaLinearizedSolver

void NonLinearWrapper::solve(){
  if(mySolver == NULL){
    throw(ErrorHandle("NonLinearWrapper", "solve", "the Solver must be set before attempting to solve"));
  }
  if(currentSolution == NULL){
    throw(ErrorHandle("NonLinearWrapper", "solve", "the current and previous Solutions should be set before attempting to solve"));
  }
  int rank;
  Partitioner * part = currentSolution->getMesh()->getPartitioner();
  bool speak = (verbose and ((part == NULL) or (part->getRank() == 0)));
  if(speak){
    std::cout << "Starting non-linear loop..." << std::endl;
  }
  for(int i = 0; i < maxIters; i++){
    linearizedSolver(mySolver);
    this->residual = residualComputer(currentSolution, previousSolution);
    if(speak){
      std::cout << "\r" << "Iteration " + std::to_string(i+1) + " at residual value " + std::to_string(this->residual) << std::flush;
    }
    if(residual < this->resTol){
      break;
    } else {
      double val = 0.0;
      for(int iVal = 0; iVal < currentSolution->getValues()->size(); iVal++){
        val = (1.0 - dampening)*(currentSolution->getValues()->at(iVal)) + dampening*(previousSolution->getValues()->at(iVal));
        currentSolution->getValues()->at(iVal) = val;
        previousSolution->getValues()->at(iVal) = val;
      }
      //std::copy(currentSolution->getValues()->begin(), currentSolution->getValues()->end(), previousSolution->getValues()->begin());
      if(part != NULL){
        part->updateSharedInformation();
      }
    }
  }
  if(speak){
    std::cout << std::endl;
    std::cout << "... finished non-linear loop." << std::endl;
  }
};//solve

};//hfox

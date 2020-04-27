#include "PetscInterface.h"

namespace hfox{

PetscInterface::PetscInterface(){
  setUp();
};//empty constructor

PetscInterface::PetscInterface(PetscOpts options){
  setUp();
  myOptions = options;
};//user options constructor

void PetscInterface::setUp(){
  PetscBool isInit;
  pErr = PetscInitialized(&isInit);
  CHKERRXX(pErr);
  if(!isInit){
    pErr = PetscInitialize(NULL, NULL, NULL, NULL);
    CHKERRXX(pErr);
    masterPetscClass = 1;
  }
};//setUp

PetscInterface::~PetscInterface(){
  destroySystem();
  if(masterPetscClass){
    pErr = PetscFinalize();
    CHKERRXX(pErr);
  }
};//destructor

void PetscInterface::destroySystem(){
  if(initialized){
    pErr = KSPDestroy(&kspSolver);
    CHKERRXX(pErr);
    pErr = MatDestroy(&M);
    CHKERRXX(pErr);
    pErr = VecDestroy(&b);
    CHKERRXX(pErr);
  }
  initialized = 0;
  configured = 0;
  allocated = 0;
  assembled = 0;
};//destroySystem

void PetscInterface::initialize(){
  if(!initialized){
    pErr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRXX(pErr);
    pErr = MatCreate(PETSC_COMM_WORLD, &M);
    CHKERRXX(pErr);
    pErr = KSPCreate(PETSC_COMM_WORLD, &kspSolver);
    CHKERRXX(pErr);
    initialized = 1;
  }
};//initialize

void PetscInterface::configure(){
  if(initialized){
    pErr = VecSetFromOptions(b);
    CHKERRXX(pErr);
    pErr = MatSetType(M, MATMPIAIJ);
    CHKERRXX(pErr);
    pErr = MatSetFromOptions(M);
    CHKERRXX(pErr);
    pErr = KSPGetPC(kspSolver, &preCond);
    CHKERRXX(pErr);
    pErr = PCSetType(preCond, myOptions.preconditionnerType);
    CHKERRXX(pErr);
    pErr = KSPSetType(kspSolver, myOptions.solverType);
    CHKERRXX(pErr);
    pErr = KSPSetTolerances(kspSolver, myOptions.rtol, PETSC_DEFAULT, PETSC_DEFAULT, myOptions.maxits);
    CHKERRXX(pErr);
    if(myOptions.verbose){
      pErr = PetscOptionsSetValue(NULL, "-ksp_monitor", NULL);
      CHKERRXX(pErr);
    }
    pErr = KSPSetFromOptions(kspSolver);
    CHKERRXX(pErr);
    configured = 1;
  } else {
    throw(ErrorHandle("PetscInterface", "configure", "the object was not initialized before being configured."));
  }
};//configure


void PetscInterface::allocate(int ndofs, const std::vector<int> * diagSparsePattern, 
        const std::vector<int> * offSparsePattern){
  if(initialized and configured){
    nDOFs = ndofs;
    pErr = MatSetSizes(M, nDOFs, nDOFs, PETSC_DETERMINE, PETSC_DETERMINE);
    CHKERRXX(pErr);
    if(diagSparsePattern == NULL){
      pErr = MatSetUp(M);
      CHKERRXX(pErr);
    } else {
      pErr = MatMPIAIJSetPreallocation(M, 0, diagSparsePattern->data(), 0, offSparsePattern->data());
      CHKERRXX(pErr);
      pErr = MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      CHKERRXX(pErr);
    }
    pErr = VecSetSizes(b, nDOFs, PETSC_DETERMINE);
    CHKERRXX(pErr);
    pErr = VecSet(b, 0.0);
    CHKERRXX(pErr);
    allocated = 1;
  } else {
    throw(ErrorHandle("PetscInterface", "allocate", "the object was not initialized and configured before being allocated."));
  }
};//allocate

void PetscInterface::addValMatrix(int i, int j, const double & val){
  if(allocated){
    pErr = MatSetValues(M, 1, &i, 1, &j, &val, ADD_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "addValMatrix", "the matrix must be allocated before adding values."));
  }
};//addValMatrix

void PetscInterface::addValsMatrix(std::vector<int> & is, std::vector<int> & js, const double * vals){
  if(allocated){
    pErr = MatSetValues(M, is.size(), is.data(), js.size(), js.data(), vals, ADD_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "addValsMatrix", "the matrix must be allocated before adding values."));
  }
};//addValsMatrix

void PetscInterface::addValRHS(int i, const double & val){
  if(allocated){
    pErr = VecSetValues(b, 1, &i, &val, ADD_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "addValRHS", "the RHS must be allocated before adding values."));
  }
};//addValRHS

void PetscInterface::addValsRHS(std::vector<int> & is, const double * vals){
  if(allocated){
    pErr = VecSetValues(b, is.size(), is.data(), vals, ADD_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "addValsRHS", "the RHS must be allocated before adding values."));
  }
};//addValsRHS

void PetscInterface::setValMatrix(int i, int j, const double & val){
  if(allocated){
    pErr = MatSetValues(M, 1, &i, 1, &j, &val, INSERT_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "setValMatrix", "the matrix must be allocated before inserting values."));
  }
};//setValMatrix


void PetscInterface::setValsMatrix(std::vector<int> & is, std::vector<int> & js, const double * vals){
  if(allocated){
    pErr = MatSetValues(M, is.size(), is.data(), js.size(), js.data(), vals, INSERT_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "setValsMatrix", "the matrix must be allocated before setting values."));
  }
};//setValsMatrix

void PetscInterface::setValRHS(int i, const double & val){
  if(allocated){
    pErr = VecSetValues(b, 1, &i, &val, INSERT_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "setValRHS", "the RHS must be allocated before setting values."));
  }
};//setValRHS

void PetscInterface::setValsRHS(std::vector<int> & is, const double * vals){
  if(allocated){
    pErr = VecSetValues(b, is.size(), is.data(), vals, INSERT_VALUES);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "setValsRHS", "the RHS must be allocated before setting values."));
  }
};//setValsRHS

void PetscInterface::zeroOutRows(std::vector<int> & is){
  if(allocated){
    pErr = MatZeroRows(M, is.size(), is.data(), 0.0, 0, 0);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "zeroOutRows", "the matrix must be allocated before zeroing out rows."));
  }
};//zeroOutRows

void PetscInterface::assemble(){
  if(allocated){
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    assembled = 1;
  } else {
    throw(ErrorHandle("PetscInterface", "assemble", "the matrix should at least be allocated before assembling."));
  }
};//assemble

void PetscInterface::solve(std::vector<double> * solution){
  if(assembled){
    KSPSetOperators(kspSolver, M, M);
    if(myOptions.verbose){
      std::cout << "PetsC resolution (rtol = " + std::to_string(myOptions.rtol) + 
        ", maxIter = " + std::to_string(myOptions.maxits) + ", nDOFs = " +
        std::to_string(nDOFs)+"): " << std::endl;
    }
    PetscInt blocksize;
    pErr = VecGetBlockSize(b, &blocksize);
    CHKERRXX(pErr);
    solution->resize(nDOFs);
    Vec sol;
    pErr = VecCreateMPIWithArray(PETSC_COMM_WORLD, blocksize, nDOFs, PETSC_DETERMINE, solution->data(), &sol);
    CHKERRXX(pErr);
    pErr = KSPSolve(kspSolver, b, sol);
    CHKERRXX(pErr);
    pErr = VecDestroy(&sol);
    CHKERRXX(pErr);
  } else {
    throw(ErrorHandle("PetscInterface", "solve", "the matrix needs to be assembled before solving."));
  }
};//solve

void PetscInterface::clearSystem(){
  if(assembled){
    MatZeroEntries(M);
    VecSet(b, 0.0);
    assemble();
  } else {
    throw(ErrorHandle("PetscInterface", "clearSystem", "the matrix needs to be assembled at least once before clearing."));
  }
};//clearSystem

}//hfox

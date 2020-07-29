#ifndef HDGSOLVEROPTS_H
#define HDGSOLVEROPTS_H

namespace hfox{

enum HDGSolverType{
  IMPLICIT,
  WEXPLICIT,
  SEXPLICIT
};//HDGSolverType

struct HDGSolverOpts{
  HDGSolverType type = IMPLICIT;
  bool verbosity = 1;
  bool doubleValuedTau = 0;
};//HDGSolverOpts

}//hfox

#endif//HDGSOLVEROPTS_H

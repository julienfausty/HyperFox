#ifndef ASSEMBLYTYPE_H
#define ASSEMBLYTYPE_H

enum UnitAssemblyType{
  None,
  Add,
  Set
};//UnitAssemblyType

struct AssemblyType{
  UnitAssemblyType matrix;
  UnitAssemblyType rhs;
};//AssemblyType

#endif//ASSEMBLYTYPE_H

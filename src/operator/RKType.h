#ifndef RKTYPE_H
#define RKTYPE_H

namespace hfox{

enum RKType{
  //Explicit
  FEuler,
  EMidpoint,
  Heun,
  Kutta3,
  Heun3,
  SSPRK3,
  RK4,
  //Implicit
  BEuler,
  IMidpoint,
  CrankNicolson,
  KS2, //Kraaijevanger and Spijker
  QZ2, //Qin and Zhang
  ALX2, //Alexander 1977
  RK43 //four stage, 3rd order
};

}//hfox

#endif//RKTYPE_H

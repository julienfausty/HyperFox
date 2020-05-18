#include "RungeKutta.h"

namespace hfox{


void RungeKutta::setUpDB(){
  butcherDB[FEuler].resize(2,2);
  butcherDB[FEuler] << 0, 0, 0, 1;
  butcherDB[EMidpoint].resize(3,3);
  butcherDB[EMidpoint] << 
    0, 0, 0, 
    0.5, 0.5, 0.0,
    0, 0, 1.0;
  butcherDB[Heun].resize(3,3);
  butcherDB[Heun] << 
    0, 0, 0, 
    1.0, 1.0, 0.0,
    0, 0.5, 0.5;
  butcherDB[Kutta3].resize(4,4);
  butcherDB[Kutta3] << 
    0, 0, 0, 0, 
    0.5, 0.5, 0.0, 0.0,
    1.0, -1.0, 2.0, 0.0,
    0.0, 1.0/6.0, 2.0/3.0, 1.0/6.0;
  butcherDB[Heun3].resize(4,4);
  butcherDB[Heun3] << 
    0, 0, 0, 0, 
    1.0/3.0, 1.0/3.0, 0.0, 0.0,
    2.0/3.0, 0.0, 2.0/3.0, 0.0,
    0.0, 1.0/4.0, 0.0, 3.0/4.0;
  butcherDB[SSPRK3].resize(4,4);
  butcherDB[SSPRK3] << 
    0, 0, 0, 0, 
    1.0, 1.0, 0.0, 0.0,
    1.0/2.0, 1.0/4.0, 1.0/4.0, 0.0,
    0.0, 1.0/6.0, 1.0/6.0, 2.0/3.0;
  butcherDB[RK4].resize(5,5);
  butcherDB[RK4] << 
    0, 0, 0, 0, 0, 
    1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0,
    1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0;
  butcherDB[BEuler].resize(2,2);
  butcherDB[BEuler] << 1, 1, 0, 1;
  butcherDB[BEuler].resize(2,2);
  butcherDB[BEuler] << 1.0/2.0, 1.0/2.0, 0, 1;
  butcherDB[CrankNicolson].resize(3,3);
  butcherDB[CrankNicolson] << 
    0.0, 0.0, 0.0,
    1.0, 1.0/2.0, 1.0/2.0,
    0.0, 1.0/2.0, 1.0/2.0;
  butcherDB[KS2].resize(3,3);
  butcherDB[KS2] << 
    1.0/2.0, 1.0/2.0, 0.0,
    3.0/2.0, -1.0/2.0, 2.0,
    0.0, -1.0/2.0, 3.0/2.0;
  butcherDB[QZ2].resize(3,3);
  butcherDB[QZ2] << 
    1.0/4.0, 1.0/4.0, 0.0,
    3.0/4.0, 1.0/2.0, 1.0/4.0,
    0.0, 1.0/2.0, 1.0/2.0;
  butcherDB[QZ2].resize(5,5);
  butcherDB[QZ2] << 
    1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0,
    2.0/3.0, 1.0/6.0, 1.0/2.0, 0.0, 0.0,
    1.0/2.0, -1.0/2.0, 1.0/2.0, 1.0/2.0, 0.0,
    1.0, 3.0/2.0, -3.0/2.0, 1.0/2.0, 1.0/2.0,
    0.0, 3.0/2.0, -3.0/2.0, 1.0/2.0, 1.0/2.0;
};//setUpDB

std::map<RKType, EMatrix> RungeKutta::butcherDB;

}//hfox

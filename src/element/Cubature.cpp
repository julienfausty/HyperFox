#include "Cubature.h"

namespace hfox{

  Cubature::Cubature(int dim, int ord, elementGeometry geom){
    if(!rulesExist){
      initializeDatabase();
      rulesExist = 1;
    }
    if(dim > Cubature::maxDim){
      ErrorHandle eh("Cubature", "Constructor", "the requested space dimension (" + std::to_string(dim)
          + ") is too large and thus not implemented in the cubature rules yet. (maxDim = " + 
          std::to_string(Cubature::maxDim) +")");
      throw(eh);
    } else if(ord > Cubature::maxOrderMap[dim]){
      ErrorHandle eh("Cubature", "Constructor", "the requested polynomial order (" + std::to_string(ord)
          + ") is too large for dimension " + std::to_string(dim) 
          + " and thus not implemented in the cubature rules yet. (maxOrder = " + 
          std::to_string(Cubature::maxOrderMap[dim]) +")");
      throw(eh);
    }
    dimension = dim;
    order = ord;
    geometry = geom;
    determineRule();
  };//Constructor

  int Cubature::getDimension() const{
    return dimension;
  };//getDimension

  int Cubature::getOrder() const{
    return order;
  };//getOrder

  elementGeometry Cubature::getGeometry() const{
    return geometry;
  };//getGeometry

  int Cubature::getNumIPs() const{
    return nIPs;
  };//getNumIPs

  const std::vector< std::vector<double> > * Cubature::getIPCoords() const{
    return &ipCoords;
  };//getIPCoords

  const std::vector<double> * Cubature::getIPWeights() const{
    return &ipWeights;
  };//getIPWeights

  void Cubature::determineRule(){
    std::tuple<int, int, elementGeometry> properties(dimension, order, geometry);
    nIPs = nIPMap[properties];
    std::tuple<std::vector< std::vector<double> >, std::vector<double> > rule;
    rule = Cubature::rulesDatabase[std::tuple<int, int, elementGeometry>(dimension, nIPs, geometry)];
    ipCoords = std::get<0>(rule);
    ipWeights = std::get<1>(rule);
  }; //determineRule

  Cubature::~Cubature(){
    ;
  };//Destructor


  //Static members
  void Cubature::initializeDatabase(){
    maxDim = 3;
    Cubature::maxOrderMap[0] = 20;
    Cubature::maxOrderMap[1] = 20;
    Cubature::maxOrderMap[2] = 20;
    Cubature::maxOrderMap[3] = 10;

    typedef std::tuple<int, int, elementGeometry> cubDataKey;
    typedef std::tuple<std::vector< std::vector<double> >, std::vector<double> > cubDataVal;
    int tempDim, tempOrder, tempNumIPs;
    std::vector< std::vector<double> > tempCoords;
    std::vector<double> tempWeights;
    std::vector<double> * tempPvec;
    int index;


    //dimension 0
    for(int j = 0; j < (maxOrderMap[0]+1); j++){
      nIPMap[cubDataKey(0, j, simplex)] = 0;
      nIPMap[cubDataKey(0, j, orthotope)] = 0;
    }
    rulesDatabase[cubDataKey(0, 0, simplex)] = cubDataVal(std::vector< std::vector<double> >(),
        std::vector<double>() );
    rulesDatabase[cubDataKey(0, 0, orthotope)] = cubDataVal(std::vector< std::vector<double> >(),
        std::vector<double>() );

    //all order 0
    std::vector<double> volumeSimplex = {0.0, 2.0, 2.0, 4.0/3.0};
    for(int i = 1; i < (maxDim + 1); i++){
      std::vector<double> pointO(i, 0.0);
      std::vector<double> weightO(1, std::pow(2,i));
      nIPMap[cubDataKey(i, 0, orthotope)] = 1;
      rulesDatabase[cubDataKey(i, 1, orthotope)] = cubDataVal(std::vector< std::vector<double> >(1, pointO), weightO);
      std::vector<double> pointS(i, -0.5);
      std::vector<double> weightS(1, volumeSimplex[i]);
      nIPMap[cubDataKey(i, 0, simplex)] = 1;
      rulesDatabase[cubDataKey(i, 1, simplex)] = cubDataVal(std::vector< std::vector<double> >(1, pointS), weightS);
    }

    //dimension 1
    tempDim = 1;
    //analytical solution for nIPs needed for polynomial order
    for(int j = 1; j < (maxOrderMap[tempDim] + 1); j++){
      nIPMap[cubDataKey(1, j, orthotope)] = std::ceil((j+1.0)/2.0);
      nIPMap[cubDataKey(1, j, simplex)] = std::ceil((j+1.0)/2.0);
    }
    //nIP 1
    tempNumIPs = 1;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, 0.0);
    tempWeights[0] = 2.0;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 2
    tempNumIPs = 2;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, -std::sqrt(1.0/3.0));
    tempCoords[1] = std::vector<double>(1, std::sqrt(1.0/3.0));
    tempWeights[0] = 1.0; tempWeights[1] = 1.0;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 3
    tempNumIPs = 3;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, 0.0);
    tempCoords[1] = std::vector<double>(1, -std::sqrt(3.0/5.0));
    tempCoords[2] = std::vector<double>(1, std::sqrt(3.0/5.0));
    tempWeights[0] = 8.0/9.0; tempWeights[1] = 5.0/9.0; tempWeights[2] = 5.0/9.0;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 4
    tempNumIPs = 4;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0)));
    tempCoords[1] = std::vector<double>(1, -std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0)));
    tempCoords[2] = std::vector<double>(1, std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0)));
    tempCoords[3] = std::vector<double>(1, -std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0)));
    tempWeights[0] = (18.0 + std::sqrt(30.0))/36.0; tempWeights[1] = tempWeights[0];
    tempWeights[2] = (18.0 - std::sqrt(30.0))/36.0; tempWeights[3] = tempWeights[2];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 5
    tempNumIPs = 5;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, 0);
    tempCoords[1] = std::vector<double>(1, (1.0/3.0)*std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0)));
    tempCoords[2] = std::vector<double>(1, -(1.0/3.0)*std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0)));
    tempCoords[3] = std::vector<double>(1, (1.0/3.0)*std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0)));
    tempCoords[4] = std::vector<double>(1, -(1.0/3.0)*std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0)));
    tempWeights[0] = 128.0/225.0;
    tempWeights[1] = (322.0 + 13.0*std::sqrt(70.0))/900.0; tempWeights[2] = tempWeights[1];
    tempWeights[3] = (322.0 - 13.0*std::sqrt(70.0))/900.0; tempWeights[4] = tempWeights[3];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 6
    tempNumIPs = 6;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, 0.661209386466264513661);
    tempCoords[1] = std::vector<double>(1, -0.661209386466264513661);
    tempCoords[2] = std::vector<double>(1, -0.2386191860831969086305);
    tempCoords[3] = std::vector<double>(1, 0.2386191860831969086305);
    tempCoords[4] = std::vector<double>(1, -0.9324695142031520278123);
    tempCoords[5] = std::vector<double>(1, 0.9324695142031520278123);
    tempWeights[0] = 0.3607615730481386075698; tempWeights[1] = tempWeights[0];
    tempWeights[2] = 0.4679139345726910473899; tempWeights[3] = tempWeights[2];
    tempWeights[4] = 0.1713244923791703450403; tempWeights[5] = tempWeights[4];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 7
    tempNumIPs = 7;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, 0.0);
    tempCoords[1] = std::vector<double>(1, 0.4058451513773972);
    tempCoords[2] = std::vector<double>(1, -0.4058451513773972);
    tempCoords[3] = std::vector<double>(1, -0.7415311855993945);
    tempCoords[4] = std::vector<double>(1, 0.7415311855993945);
    tempCoords[5] = std::vector<double>(1, -0.9491079123427585);
    tempCoords[6] = std::vector<double>(1, 0.9491079123427585);
    tempWeights[0] = 0.4179591836734694;
    tempWeights[1] = 0.3818300505051189; tempWeights[2] = tempWeights[1];
    tempWeights[3] = 0.2797053914892766; tempWeights[4] = tempWeights[3];
    tempWeights[5] = 0.1294849661688697; tempWeights[6] = tempWeights[5];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 8
    tempNumIPs = 8;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, -0.1834346424956498);
    tempCoords[1] = std::vector<double>(1, 0.1834346424956498);
    tempCoords[2] = std::vector<double>(1, -0.5255324099163290);
    tempCoords[3] = std::vector<double>(1, 0.5255324099163290);
    tempCoords[4] = std::vector<double>(1, -0.7966664774136267);
    tempCoords[5] = std::vector<double>(1, 0.7966664774136267);
    tempCoords[6] = std::vector<double>(1, -0.9602898564975363);
    tempCoords[7] = std::vector<double>(1, 0.9602898564975363);
    tempWeights[0] = 0.3626837833783620; tempWeights[1] = tempWeights[0];
    tempWeights[2] = 0.3137066458778873; tempWeights[3] = tempWeights[2];
    tempWeights[4] = 0.2223810344533745; tempWeights[5] = tempWeights[4];
    tempWeights[6] = 0.1012285362903763; tempWeights[7] = tempWeights[6];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 9
    tempNumIPs = 9;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, 0.0);
    tempCoords[0] = std::vector<double>(1, -0.8360311073266358);
    tempCoords[1] = std::vector<double>(1, 0.8360311073266358);
    tempCoords[2] = std::vector<double>(1, -0.9681602395076261);
    tempCoords[3] = std::vector<double>(1, 0.9681602395076261);
    tempCoords[4] = std::vector<double>(1, -0.3242534234038089);
    tempCoords[5] = std::vector<double>(1, 0.3242534234038089);
    tempCoords[6] = std::vector<double>(1, -0.6133714327005904);
    tempCoords[7] = std::vector<double>(1, 0.6133714327005904);
    tempWeights[0] = 0.3302393550012598;
    tempWeights[1] = 0.1806481606948574; tempWeights[2] = tempWeights[1];
    tempWeights[3] = 0.0812743883615744; tempWeights[4] = tempWeights[3];
    tempWeights[5] = 0.3123470770400029; tempWeights[6] = tempWeights[5];
    tempWeights[7] = 0.2606106964029354; tempWeights[8] = tempWeights[7];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 10
    tempNumIPs = 10;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, -0.1488743389816312);
    tempCoords[1] = std::vector<double>(1, 0.1488743389816312);
    tempCoords[2] = std::vector<double>(1, -0.4333953941292472);
    tempCoords[3] = std::vector<double>(1, 0.4333953941292472);
    tempCoords[4] = std::vector<double>(1, -0.6794095682990244);
    tempCoords[5] = std::vector<double>(1, 0.6794095682990244);
    tempCoords[6] = std::vector<double>(1, -0.8650633666889845);
    tempCoords[7] = std::vector<double>(1, 0.8650633666889845);
    tempCoords[8] = std::vector<double>(1, -0.9739065285171717);
    tempCoords[9] = std::vector<double>(1, 0.9739065285171717);
    tempWeights[0] = 0.2955242247147529; tempWeights[1] = tempWeights[0];
    tempWeights[2] = 0.2692667193099963; tempWeights[3] = tempWeights[2];
    tempWeights[4] = 0.2190863625159820; tempWeights[5] = tempWeights[4];
    tempWeights[6] = 0.1494513491505806; tempWeights[7] = tempWeights[6];
    tempWeights[8] = 0.0666713443086881; tempWeights[9] = tempWeights[8];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //nIP 11
    tempNumIPs = 11;
    tempCoords.resize(tempNumIPs);
    tempWeights.resize(tempNumIPs);
    tempCoords[0] = std::vector<double>(1, 0.0);
    tempCoords[1] = std::vector<double>(1, -0.2695431559523450);
    tempCoords[2] = std::vector<double>(1, 0.2695431559523450);
    tempCoords[3] = std::vector<double>(1, -0.5190961292068118);
    tempCoords[4] = std::vector<double>(1, 0.5190961292068118);
    tempCoords[5] = std::vector<double>(1, -0.7301520055740494);
    tempCoords[6] = std::vector<double>(1, 0.7301520055740494);
    tempCoords[7] = std::vector<double>(1, -0.8870625997680953);
    tempCoords[8] = std::vector<double>(1, 0.8870625997680953);
    tempCoords[9] = std::vector<double>(1, -0.9782286581460570);
    tempCoords[10] = std::vector<double>(1, 0.9782286581460570);
    tempWeights[0] = 0.2729250867779006;
    tempWeights[1] = 0.2628045445102467; tempWeights[2] = tempWeights[1];
    tempWeights[3] = 0.2331937645919905; tempWeights[4] = tempWeights[3];
    tempWeights[5] = 0.1862902109277343; tempWeights[6] = tempWeights[5];
    tempWeights[7] = 0.1255803694649046; tempWeights[8] = tempWeights[7];
    tempWeights[9] = 0.0556685671161737; tempWeights[10] = tempWeights[9];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    // dimension 2
    tempDim = 2;
    // order 1
    tempNumIPs = 1;
    nIPMap[cubDataKey(2, 1, orthotope)] = 1;
    nIPMap[cubDataKey(2, 1, simplex)] = 1;
    // nIP 1
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempCoords[0][0] = 0.0; tempCoords[0][1] = 0.0;
    tempWeights[0] = 4.0;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    tempCoords[0][0] = -1.0/3.0; tempCoords[0][1] = -1.0/3.0;
    tempWeights[0] = 2.0;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);

    //order 2 simplex
    tempNumIPs = 3;
    nIPMap[cubDataKey(2, 2, simplex)] = tempNumIPs;
    // nIP 3
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      (tempCoords[i]).resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempCoords[0][0] = -2.0/3.0; tempCoords[0][1] = 1.0/3.0;
    tempCoords[1][0] = 1.0/3.0; tempCoords[1][1] = -2.0/3.0;
    tempCoords[2][0] = -2.0/3.0; tempCoords[2][1] = -2.0/3.0;
    tempWeights[0] = 2.0/3.0; tempWeights[1] = tempWeights[0]; tempWeights[2] = tempWeights[0];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);

    //order 2 and 3 orthotope
    tempNumIPs = 4;
    nIPMap[cubDataKey(2, 2, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 3, orthotope)] = tempNumIPs;
    // nIP 4
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempCoords[0][0] = 0.57735026918962576450914878050195745565; tempCoords[0][1] = tempCoords[0][0];
    tempCoords[1][0] = tempCoords[0][0]; tempCoords[1][1] = -tempCoords[0][0];
    tempCoords[2][0] = -tempCoords[0][0]; tempCoords[2][1] = tempCoords[0][0];
    tempCoords[3][0] = -tempCoords[0][0]; tempCoords[3][1] = -tempCoords[0][0];
    tempWeights[0] = 1.0; tempWeights[1] = tempWeights[0]; 
    tempWeights[2] = tempWeights[0]; tempWeights[3] = tempWeights[0];
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);

    //order 3 and 4 simplex
    tempNumIPs = 6;
    nIPMap[cubDataKey(2, 3, simplex)] = tempNumIPs;
    nIPMap[cubDataKey(2, 4, simplex)] = tempNumIPs;
    // nIP 6
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
{-0.1081030181680702273633414922338960232,   -0.7837939636638595452733170155322079536,   0.44676317935602293139001401686624560874,
 -0.7837939636638595452733170155322079536,   -0.1081030181680702273633414922338960232,   0.44676317935602293139001401686624560874,
 -0.1081030181680702273633414922338960232,   -0.1081030181680702273633414922338960232,  0.44676317935602293139001401686624560874,
-0.81684757298045851308085707319559698429,   0.63369514596091702616171414639119396858,   0.21990348731064373527665264980042105793,
 0.63369514596091702616171414639119396858,  -0.81684757298045851308085707319559698429,   0.21990348731064373527665264980042105793,
-0.81684757298045851308085707319559698429,  -0.81684757298045851308085707319559698429,   0.21990348731064373527665264980042105793});
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);

    //order 4 and 5 orthotope
    tempNumIPs = 8;
    nIPMap[cubDataKey(2, 4, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 5, orthotope)] = tempNumIPs;
    // nIP 8
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0.68313005106397322554806924536807013272,                                          0,    0.8163265306122448979591836734693877551,
                                        0,   0.68313005106397322554806924536807013272,    0.8163265306122448979591836734693877551,
-0.68313005106397322554806924536807013272,                                          0,    0.8163265306122448979591836734693877551,
                                        0,  -0.68313005106397322554806924536807013272,    0.8163265306122448979591836734693877551,
  0.8819171036881968635005385845464201419,    0.8819171036881968635005385845464201419,    0.1836734693877551020408163265306122449,
  0.8819171036881968635005385845464201419,   -0.8819171036881968635005385845464201419,    0.1836734693877551020408163265306122449,
 -0.8819171036881968635005385845464201419,    0.8819171036881968635005385845464201419,    0.1836734693877551020408163265306122449,
 -0.8819171036881968635005385845464201419,   -0.8819171036881968635005385845464201419,    0.1836734693877551020408163265306122449}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 5 simplex
    tempNumIPs = 7;
    nIPMap[cubDataKey(2, 5, simplex)] = tempNumIPs;
    // nIP 7
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        { -0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,                                        0.45,
 -0.79742698535308732239802527616975234389,    0.59485397070617464479605055233950468778,    0.25187836108965430519136789100036266732,
  0.59485397070617464479605055233950468778,   -0.79742698535308732239802527616975234389,    0.25187836108965430519136789100036266732,
 -0.79742698535308732239802527616975234389,   -0.79742698535308732239802527616975234389,    0.25187836108965430519136789100036266732,
-0.059715871789769820459117580973104798968,   -0.88056825642046035908176483805379040206,    0.26478830557701236147529877566630399935,
 -0.88056825642046035908176483805379040206,  -0.059715871789769820459117580973104798968,    0.26478830557701236147529877566630399935,
-0.059715871789769820459117580973104798968,  -0.059715871789769820459117580973104798968,    0.26478830557701236147529877566630399935}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 6 and 7 orthotope
    tempNumIPs = 12;
    nIPMap[cubDataKey(2, 6, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 7, orthotope)] = tempNumIPs;
    // nIP 12
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0.92582009977255146156656677658399952253,                                          0,   0.24197530864197530864197530864197530864,
                                        0,   0.92582009977255146156656677658399952253,   0.24197530864197530864197530864197530864,
-0.92582009977255146156656677658399952253,                                          0,   0.24197530864197530864197530864197530864,
                                        0,  -0.92582009977255146156656677658399952253,   0.24197530864197530864197530864197530864,
  0.8059797829185987437078561813507442463,    0.8059797829185987437078561813507442463,   0.23743177469063023421810525931129352533,
  0.8059797829185987437078561813507442463,   -0.8059797829185987437078561813507442463,   0.23743177469063023421810525931129352533,
 -0.8059797829185987437078561813507442463,    0.8059797829185987437078561813507442463,   0.23743177469063023421810525931129352533,
 -0.8059797829185987437078561813507442463,   -0.8059797829185987437078561813507442463,   0.23743177469063023421810525931129352533,
  0.3805544332083156563791063590863941355,    0.3805544332083156563791063590863941355,   0.52059291666739445713991943204673116603,
  0.3805544332083156563791063590863941355,   -0.3805544332083156563791063590863941355,   0.52059291666739445713991943204673116603,
 -0.3805544332083156563791063590863941355,    0.3805544332083156563791063590863941355,   0.52059291666739445713991943204673116603,
 -0.3805544332083156563791063590863941355,   -0.3805544332083156563791063590863941355,   0.52059291666739445713991943204673116603}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 6 simplex
    tempNumIPs = 12;
    nIPMap[cubDataKey(2, 6, simplex)] = tempNumIPs;
    // nIP 12
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.87382197101699554331933679425836168532 ,   0.74764394203399108663867358851672337064,    0.10168981274041363384187361821373796809,
  0.74764394203399108663867358851672337064,   -0.87382197101699554331933679425836168532,    0.10168981274041363384187361821373796809,
 -0.87382197101699554331933679425836168532,   -0.87382197101699554331933679425836168532,    0.10168981274041363384187361821373796809,
 -0.50142650965817915741672289378596184782,  0.0028530193163583148334457875719236956481,    0.23357255145275873205057922277115888265,
0.0028530193163583148334457875719236956481,   -0.50142650965817915741672289378596184782,    0.23357255145275873205057922277115888265,
 -0.50142650965817915741672289378596184782,   -0.50142650965817915741672289378596184782,    0.23357255145275873205057922277115888265,
 -0.89370990031036610529350065673720370601,    0.27300499824279729446028518882409939961,    0.16570215123674715038710691284088490796,
  0.27300499824279729446028518882409939961,   -0.89370990031036610529350065673720370601,    0.16570215123674715038710691284088490796,
 -0.37929509793243118916678453208689569359,    0.27300499824279729446028518882409939961,    0.16570215123674715038710691284088490796,
  0.27300499824279729446028518882409939961,   -0.37929509793243118916678453208689569359,    0.16570215123674715038710691284088490796,
 -0.37929509793243118916678453208689569359,   -0.89370990031036610529350065673720370601,    0.16570215123674715038710691284088490796,
 -0.89370990031036610529350065673720370601,   -0.37929509793243118916678453208689569359,    0.16570215123674715038710691284088490796}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 7 simplex
    tempNumIPs = 15;
    nIPMap[cubDataKey(2, 7, simplex)] = tempNumIPs;
    // nIP 15
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.93253870289082430257005654739836752385,     0.8650774057816486051401130947967350477,   0.033090100221584262071955896945834891238,
   0.8650774057816486051401130947967350477,   -0.93253870289082430257005654739836752385,   0.033090100221584262071955896945834891238,
 -0.93253870289082430257005654739836752385,   -0.93253870289082430257005654739836752385,   0.033090100221584262071955896945834891238,
 -0.51684523480919288209962646032443600082,   0.033690469618385764199252920648872001638,    0.25588834246031114556580247036929262133,
 0.033690469618385764199252920648872001638,   -0.51684523480919288209962646032443600082,    0.25588834246031114556580247036929262133,
 -0.51684523480919288209962646032443600082,   -0.51684523480919288209962646032443600082,    0.25588834246031114556580247036929262133,
-0.051380614990563531580838528101628435036,   -0.89723877001887293683832294379674312993,    0.15417329237197213566964304166748277293,
 -0.89723877001887293683832294379674312993,  -0.051380614990563531580838528101628435036,    0.15417329237197213566964304166748277293,
-0.051380614990563531580838528101628435036,  -0.051380614990563531580838528101628435036,    0.15417329237197213566964304166748277293,
 -0.90592671069480953331718004928623004897,    0.50856008110010635471247864925623992738,    0.11175746580639956167963262884202819059,
  0.50856008110010635471247864925623992738,   -0.90592671069480953331718004928623004897,    0.11175746580639956167963262884202819059,
 -0.60263337040529682139529859997000987842,    0.50856008110010635471247864925623992738,    0.11175746580639956167963262884202819059,
  0.50856008110010635471247864925623992738,   -0.60263337040529682139529859997000987842,    0.11175746580639956167963262884202819059,
 -0.60263337040529682139529859997000987842,   -0.90592671069480953331718004928623004897,    0.11175746580639956167963262884202819059,
 -0.90592671069480953331718004928623004897,   -0.60263337040529682139529859997000987842,    0.11175746580639956167963262884202819059}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 8 simplex
    tempNumIPs = 16;
    nIPMap[cubDataKey(2, 8, simplex)] = tempNumIPs;
    // nIP 16
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,     0.2886312153555743365021822209781292496,
-0.081414823414553687942368971011661355879,   -0.83717035317089262411526205797667728824,     0.1901832685345692495877922087771686332,
 -0.83717035317089262411526205797667728824,  -0.081414823414553687942368971011661355879,     0.1901832685345692495877922087771686332,
-0.081414823414553687942368971011661355879,  -0.081414823414553687942368971011661355879,     0.1901832685345692495877922087771686332,
 -0.65886138449647958675541299701707099796,    0.31772276899295917351082599403414199593,    0.20643474106943650056358310058425806003,
  0.31772276899295917351082599403414199593,   -0.65886138449647958675541299701707099796,    0.20643474106943650056358310058425806003,
 -0.65886138449647958675541299701707099796,   -0.65886138449647958675541299701707099796,    0.20643474106943650056358310058425806003,
 -0.89890554336593804908315289880680210631,    0.79781108673187609816630579761360421262,   0.064916995246396160621851856683561193593,
  0.79781108673187609816630579761360421262,   -0.89890554336593804908315289880680210631,   0.064916995246396160621851856683561193593,
 -0.89890554336593804908315289880680210631,   -0.89890554336593804908315289880680210631,   0.064916995246396160621851856683561193593,
 -0.98321044518008478932557233092141110162,    0.45698478591080856248200075835212392604,    0.05446062834886998852968938014781784832,
  0.45698478591080856248200075835212392604,   -0.98321044518008478932557233092141110162,    0.05446062834886998852968938014781784832,
 -0.47377434073072377315642842743071282442,    0.45698478591080856248200075835212392604,    0.05446062834886998852968938014781784832,
  0.45698478591080856248200075835212392604,   -0.47377434073072377315642842743071282442,    0.05446062834886998852968938014781784832,
 -0.47377434073072377315642842743071282442,   -0.98321044518008478932557233092141110162,    0.05446062834886998852968938014781784832,
 -0.98321044518008478932557233092141110162,   -0.47377434073072377315642842743071282442,    0.05446062834886998852968938014781784832}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 9 simplex
    tempNumIPs = 19;
    nIPMap[cubDataKey(2, 9, simplex)] = tempNumIPs;
    // nIP 19
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        { -0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,    0.19427159256559766763848396501457725948,
 -0.12582081701412672546013927112929005464,   -0.74835836597174654907972145774141989072,    0.15565508200954855863347871259880791224,
 -0.74835836597174654907972145774141989072,   -0.12582081701412672546013927112929005464,    0.15565508200954855863347871259880791224,
 -0.12582081701412672546013927112929005464,   -0.12582081701412672546013927112929005464,    0.15565508200954855863347871259880791224,
 -0.62359292876193453951807743906532886378,    0.24718585752386907903615487813065772755,     0.1592954778544205060657835485280905465,
  0.24718585752386907903615487813065772755,   -0.62359292876193453951807743906532886378,     0.1592954778544205060657835485280905465,
 -0.62359292876193453951807743906532886378,   -0.62359292876193453951807743906532886378,     0.1592954778544205060657835485280905465,
-0.020634961602524744432586150327614401262,   -0.95873007679495051113482769934477119748,   0.062669400454278141073709662574418630986,
 -0.95873007679495051113482769934477119748,  -0.020634961602524744432586150327614401262,   0.062669400454278141073709662574418630986,
-0.020634961602524744432586150327614401262,  -0.020634961602524744432586150327614401262,   0.062669400454278141073709662574418630986,
    -0.91054097321109458026978682006744727,       0.82108194642218916053957364013489454,   0.051155351317396062523357597117999647955,
     0.82108194642218916053957364013489454,      -0.91054097321109458026978682006744727,   0.051155351317396062523357597117999647955,
    -0.91054097321109458026978682006744727,      -0.91054097321109458026978682006744727,   0.051155351317396062523357597117999647955,
 -0.92632317589052743273036480243322979538,    0.48239719756899604138015974704684765487,   0.086567078754578754578754578754578754579,
  0.48239719756899604138015974704684765487,   -0.92632317589052743273036480243322979538,   0.086567078754578754578754578754578754579,
 -0.55607402167846860864979494461361785949,    0.48239719756899604138015974704684765487,   0.086567078754578754578754578754578754579,
  0.48239719756899604138015974704684765487,   -0.55607402167846860864979494461361785949,   0.086567078754578754578754578754578754579,
 -0.55607402167846860864979494461361785949,   -0.92632317589052743273036480243322979538,   0.086567078754578754578754578754578754579,
 -0.92632317589052743273036480243322979538,   -0.55607402167846860864979494461361785949,   0.086567078754578754578754578754578754579}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 8 and 9 orthotope
    tempNumIPs = 20;
    nIPMap[cubDataKey(2, 8, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 9, orthotope)] = tempNumIPs;
    // nIP 20
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        { 0.48892685697436906957409306022929339497,                                          0,   0.45416396068674902744178724560528503088,
                                        0,   0.48892685697436906957409306022929339497,   0.45416396068674902744178724560528503088,
-0.48892685697436906957409306022929339497,                                          0,   0.45416396068674902744178724560528503088,
                                        0,  -0.48892685697436906957409306022929339497,   0.45416396068674902744178724560528503088,
 0.93965525809683770585757175046316803289,   0.93965525809683770585757175046316803289,  0.042731231865775762580130721539841651755,
 0.93965525809683770585757175046316803289,  -0.93965525809683770585757175046316803289,  0.042731231865775762580130721539841651755,
-0.93965525809683770585757175046316803289,   0.93965525809683770585757175046316803289,  0.042731231865775762580130721539841651755,
-0.93965525809683770585757175046316803289,  -0.93965525809683770585757175046316803289,  0.042731231865775762580130721539841651755,
 0.69088055048634387281488830813840885584,   0.69088055048634387281488830813840885584,   0.21420036092686163194242563430597811641,
 0.69088055048634387281488830813840885584,  -0.69088055048634387281488830813840885584,   0.21420036092686163194242563430597811641,
-0.69088055048634387281488830813840885584,   0.69088055048634387281488830813840885584,   0.21420036092686163194242563430597811641,
-0.69088055048634387281488830813840885584,  -0.69088055048634387281488830813840885584,   0.21420036092686163194242563430597811641,
  0.9186204410567222596554701791427794163,    0.3448720253644035761712304181246973037,   0.14445222326030678901782819927444760048,
  0.3448720253644035761712304181246973037,    0.9186204410567222596554701791427794163,   0.14445222326030678901782819927444760048,
  0.9186204410567222596554701791427794163,   -0.3448720253644035761712304181246973037,   0.14445222326030678901782819927444760048,
 -0.3448720253644035761712304181246973037,    0.9186204410567222596554701791427794163,   0.14445222326030678901782819927444760048,
 -0.9186204410567222596554701791427794163,    0.3448720253644035761712304181246973037,   0.14445222326030678901782819927444760048,
  0.3448720253644035761712304181246973037,   -0.9186204410567222596554701791427794163,   0.14445222326030678901782819927444760048,
 -0.9186204410567222596554701791427794163,   -0.3448720253644035761712304181246973037,   0.14445222326030678901782819927444760048,
 -0.3448720253644035761712304181246973037,   -0.9186204410567222596554701791427794163,   0.14445222326030678901782819927444760048}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 10 simplex
    tempNumIPs = 25;
    nIPMap[cubDataKey(2, 10, simplex)] = tempNumIPs;
    // nIP 25
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,  -0.33333333333333333333333333333333333333,   0.16348665829257193285623736996835522439,
-0.93588925356611297413803082132702052207,   0.87177850713222594827606164265404104414,  0.026705937626299132551145956798137306994,
 0.87177850713222594827606164265404104414,  -0.93588925356611297413803082132702052207,  0.026705937626299132551145956798137306994,
-0.93588925356611297413803082132702052207,  -0.93588925356611297413803082132702052207,  0.026705937626299132551145956798137306994,
-0.71567779788687122981567579361808335656,   0.43135559577374245963135158723616671312,    0.0919159272094894560275758192650956328,
 0.43135559577374245963135158723616671312,  -0.71567779788687122981567579361808335656,    0.0919159272094894560275758192650956328,
-0.71567779788687122981567579361808335656,  -0.71567779788687122981567579361808335656,    0.0919159272094894560275758192650956328,
-0.35637400942232915754980487802790257253,  0.060108237854688056554191347891388233782,   0.12780981279284809086579746752530664373,
0.060108237854688056554191347891388233782,  -0.35637400942232915754980487802790257253,   0.12780981279284809086579746752530664373,
-0.70373422843235889900438646986348566125,  0.060108237854688056554191347891388233782,   0.12780981279284809086579746752530664373,
0.060108237854688056554191347891388233782,  -0.70373422843235889900438646986348566125,   0.12780981279284809086579746752530664373,
-0.70373422843235889900438646986348566125,  -0.35637400942232915754980487802790257253,   0.12780981279284809086579746752530664373,
-0.35637400942232915754980487802790257253,  -0.70373422843235889900438646986348566125,   0.12780981279284809086579746752530664373,
-0.94076022102254046473232746114791439895,   0.20246665736691849090948578691737561875,  0.068369296325918857257383168082684580517,
 0.20246665736691849090948578691737561875,  -0.94076022102254046473232746114791439895,  0.068369296325918857257383168082684580517,
 -0.2617064363443780261771583257694612198,   0.20246665736691849090948578691737561875,  0.068369296325918857257383168082684580517,
 0.20246665736691849090948578691737561875,   -0.2617064363443780261771583257694612198,  0.068369296325918857257383168082684580517,
 -0.2617064363443780261771583257694612198,  -0.94076022102254046473232746114791439895,  0.068369296325918857257383168082684580517,
-0.94076022102254046473232746114791439895,   -0.2617064363443780261771583257694612198,  0.068369296325918857257383168082684580517,
-0.94326466932012312149912848884373958881,   0.61586120184575813015989980576348821814,  0.050595515414576768778085581365666435125,
 0.61586120184575813015989980576348821814,  -0.94326466932012312149912848884373958881,  0.050595515414576768778085581365666435125,
-0.67259653252563500866077131691974862932,   0.61586120184575813015989980576348821814,  0.050595515414576768778085581365666435125,
 0.61586120184575813015989980576348821814,  -0.67259653252563500866077131691974862932,  0.050595515414576768778085581365666435125,
-0.67259653252563500866077131691974862932,  -0.94326466932012312149912848884373958881,  0.050595515414576768778085581365666435125,
-0.94326466932012312149912848884373958881,  -0.67259653252563500866077131691974862932,  0.050595515414576768778085581365666435125}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 10 and 11 orthotope
    tempNumIPs = 28;
    nIPMap[cubDataKey(2, 10, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 11, orthotope)] = tempNumIPs;
    // nIP 28
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0.71461782966460591762942382404176367266,                                          0,   0.21740043986871200551566672515052291339,
                                        0,   0.71461782966460591762942382404176367266,   0.21740043986871200551566672515052291339,
-0.71461782966460591762942382404176367266,                                          0,   0.21740043986871200551566672515052291339,
                                        0,  -0.71461782966460591762942382404176367266,   0.21740043986871200551566672515052291339,
 0.27365721017145961383557026431906649027,   0.27365721017145961383557026431906649027,   0.27727410298385108795028691796804462632,
 0.27365721017145961383557026431906649027,  -0.27365721017145961383557026431906649027,   0.27727410298385108795028691796804462632,
-0.27365721017145961383557026431906649027,   0.27365721017145961383557026431906649027,   0.27727410298385108795028691796804462632,
-0.27365721017145961383557026431906649027,  -0.27365721017145961383557026431906649027,   0.27727410298385108795028691796804462632,
 0.63660393221230104405852833198974404207,   0.63660393221230104405852833198974404207,   0.21393363787824810450280500660206134837,
 0.63660393221230104405852833198974404207,  -0.63660393221230104405852833198974404207,   0.21393363787824810450280500660206134837,
-0.63660393221230104405852833198974404207,   0.63660393221230104405852833198974404207,   0.21393363787824810450280500660206134837,
-0.63660393221230104405852833198974404207,  -0.63660393221230104405852833198974404207,   0.21393363787824810450280500660206134837,
 0.95163038878403345881049829363696371918,   0.81556543368963841306389865463928759233,  0.044074569114983092054304271154112061776,
 0.81556543368963841306389865463928759233,   0.95163038878403345881049829363696371918,  0.044074569114983092054304271154112061776,
 0.95163038878403345881049829363696371918,  -0.81556543368963841306389865463928759233,  0.044074569114983092054304271154112061776,
-0.81556543368963841306389865463928759233,   0.95163038878403345881049829363696371918,  0.044074569114983092054304271154112061776,
-0.95163038878403345881049829363696371918,   0.81556543368963841306389865463928759233,  0.044074569114983092054304271154112061776,
 0.81556543368963841306389865463928759233,  -0.95163038878403345881049829363696371918,  0.044074569114983092054304271154112061776,
-0.95163038878403345881049829363696371918,  -0.81556543368963841306389865463928759233,  0.044074569114983092054304271154112061776,
-0.81556543368963841306389865463928759233,  -0.95163038878403345881049829363696371918,  0.044074569114983092054304271154112061776,
 0.34620720004764544118747320724330043979,   0.93556787148759108135480212161830515337,   0.10162134051961130896131640398557349419,
 0.93556787148759108135480212161830515337,   0.34620720004764544118747320724330043979,   0.10162134051961130896131640398557349419,
 0.34620720004764544118747320724330043979,  -0.93556787148759108135480212161830515337,   0.10162134051961130896131640398557349419,
-0.93556787148759108135480212161830515337,   0.34620720004764544118747320724330043979,   0.10162134051961130896131640398557349419,
-0.34620720004764544118747320724330043979,   0.93556787148759108135480212161830515337,   0.10162134051961130896131640398557349419,
 0.93556787148759108135480212161830515337,  -0.34620720004764544118747320724330043979,   0.10162134051961130896131640398557349419,
-0.34620720004764544118747320724330043979,  -0.93556787148759108135480212161830515337,   0.10162134051961130896131640398557349419,
-0.93556787148759108135480212161830515337,  -0.34620720004764544118747320724330043979,   0.10162134051961130896131640398557349419}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 11 simplex
    tempNumIPs = 28;
    nIPMap[cubDataKey(2, 11, simplex)] = tempNumIPs;
    // nIP 28
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,    -0.33333333333333333333333333333333333333,     0.17152235946444843427639825073065398895,
  -0.94302916477125618104036726638981901865,      0.8860583295425123620807345327796380373,    0.020863741025789391746056984197324633703,
    0.8860583295425123620807345327796380373,    -0.94302916477125618104036726638981901865,    0.020863741025789391746056984197324633703,
  -0.94302916477125618104036726638981901865,    -0.94302916477125618104036726638981901865,    0.020863741025789391746056984197324633703,
  -0.57956008659364348598387546330841477973,     0.15912017318728697196775092661682955946,     0.14103136822343315668414771834585673963,
   0.15912017318728697196775092661682955946,    -0.57956008659364348598387546330841477973,     0.14103136822343315668414771834585673963,
  -0.57956008659364348598387546330841477973,    -0.57956008659364348598387546330841477973,     0.14103136822343315668414771834585673963,
  -0.79472903457550713910176986023011167594,     0.58945806915101427820353972046022335189,    0.077261518474038644953267609178802298293,
   0.58945806915101427820353972046022335189,    -0.79472903457550713910176986023011167594,    0.077261518474038644953267609178802298293,
  -0.79472903457550713910176986023011167594,    -0.79472903457550713910176986023011167594,    0.077261518474038644953267609178802298293,
-0.0082161980682181738701197619927680147278,    -0.98356760386356365225976047601446397054,     0.03321254610917073915275167855011212173,
  -0.98356760386356365225976047601446397054,  -0.0082161980682181738701197619927680147278,     0.03321254610917073915275167855011212173,
-0.0082161980682181738701197619927680147278,  -0.0082161980682181738701197619927680147278,     0.03321254610917073915275167855011212173,
  -0.12306814647129561795332827431496109166,    -0.75386370705740876409334345137007781667,     0.13463230815893660239567441124657557177,
  -0.75386370705740876409334345137007781667,    -0.12306814647129561795332827431496109166,     0.13463230815893660239567441124657557177,
  -0.12306814647129561795332827431496109166,    -0.12306814647129561795332827431496109166,     0.13463230815893660239567441124657557177,
  -0.70135042269583522760854842424565295474,     0.68669956732370632439798287279764986108,    0.020580579145906554986161254967959059438,
   0.68669956732370632439798287279764986108,    -0.70135042269583522760854842424565295474,    0.020580579145906554986161254967959059438,
  -0.98534914462787109678943444855199690633,     0.68669956732370632439798287279764986108,    0.020580579145906554986161254967959059438,
   0.68669956732370632439798287279764986108,    -0.98534914462787109678943444855199690633,    0.020580579145906554986161254967959059438,
  -0.98534914462787109678943444855199690633,    -0.70135042269583522760854842424565295474,    0.020580579145906554986161254967959059438,
  -0.70135042269583522760854842424565295474,    -0.98534914462787109678943444855199690633,    0.020580579145906554986161254967959059438,
  -0.90797899966914008822008044741026548976,      0.3288167483937283948070953526650315723,    0.080664953281001105168489835817596259842,
    0.3288167483937283948070953526650315723,    -0.90797899966914008822008044741026548976,    0.080664953281001105168489835817596259842,
  -0.42083774872458830658701490525476608254,      0.3288167483937283948070953526650315723,    0.080664953281001105168489835817596259842,
    0.3288167483937283948070953526650315723,    -0.42083774872458830658701490525476608254,    0.080664953281001105168489835817596259842,
  -0.42083774872458830658701490525476608254,    -0.90797899966914008822008044741026548976,    0.080664953281001105168489835817596259842,
  -0.90797899966914008822008044741026548976,    -0.42083774872458830658701490525476608254,    0.080664953281001105168489835817596259842}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 12 simplex
    tempNumIPs = 33;
    nIPMap[cubDataKey(2, 12, simplex)] = tempNumIPs;
    // nIP 33
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.023592498108916896443983752262370486151,    -0.9528150037821662071120324954752590277,   0.048533676162904066301435945151216167798,
  -0.9528150037821662071120324954752590277,  -0.023592498108916896443983752262370486151,   0.048533676162904066301435945151216167798,
-0.023592498108916896443983752262370486151,  -0.023592498108916896443983752262370486151,   0.048533676162904066301435945151216167798,
 -0.78148434468129141883250537599138250938,    0.56296868936258283766501075198276501876,   0.056972104137755089999388947815030892501,
  0.56296868936258283766501075198276501876,   -0.78148434468129141883250537599138250938,   0.056972104137755089999388947815030892501,
 -0.78148434468129141883250537599138250938,   -0.78148434468129141883250537599138250938,   0.056972104137755089999388947815030892501,
 -0.45707498597014783024400810395490980293,   -0.08585002805970433951198379209018039414,    0.12508242639180552093864756289103967826,
 -0.08585002805970433951198379209018039414,   -0.45707498597014783024400810395490980293,    0.12508242639180552093864756289103967826,
 -0.45707498597014783024400810395490980293,   -0.45707498597014783024400810395490980293,    0.12508242639180552093864756289103967826,
  -0.9507072731273288104665376055713308436,    0.90141454625465762093307521114266168719,   0.015863285019947276918629575811719468968,
  0.90141454625465762093307521114266168719,    -0.9507072731273288104665376055713308436,   0.015863285019947276918629575811719468968,
  -0.9507072731273288104665376055713308436,    -0.9507072731273288104665376055713308436,   0.015863285019947276918629575811719468968,
 -0.11977670268281377797309259385418361789,   -0.76044659463437244405381481229163276421,   0.099836669856121884238103444121988601916,
 -0.76044659463437244405381481229163276421,   -0.11977670268281377797309259385418361789,   0.099836669856121884238103444121988601916,
 -0.11977670268281377797309259385418361789,   -0.11977670268281377797309259385418361789,   0.099836669856121884238103444121988601916,
 -0.41668864052331807893217237613211876563,     0.3706203278127837999688893722742605014,   0.043567170077215115865262872678583643728,
   0.3706203278127837999688893722742605014,   -0.41668864052331807893217237613211876563,   0.043567170077215115865262872678583643728,
 -0.95393168728946572103671699614214173577,     0.3706203278127837999688893722742605014,   0.043567170077215115865262872678583643728,
   0.3706203278127837999688893722742605014,   -0.95393168728946572103671699614214173577,   0.043567170077215115865262872678583643728,
 -0.95393168728946572103671699614214173577,   -0.41668864052331807893217237613211876563,   0.043567170077215115865262872678583643728,
 -0.41668864052331807893217237613211876563,   -0.95393168728946572103671699614214173577,   0.043567170077215115865262872678583643728,
  -0.7674079606441468267385037317005977222,    0.25649950336711213367579227673445839675,    0.08645472731882842109826654183371845927,
  0.25649950336711213367579227673445839675,    -0.7674079606441468267385037317005977222,    0.08645472731882842109826654183371845927,
 -0.48909154272296530693728854503386067454,    0.25649950336711213367579227673445839675,    0.08645472731882842109826654183371845927,
  0.25649950336711213367579227673445839675,   -0.48909154272296530693728854503386067454,    0.08645472731882842109826654183371845927,
 -0.48909154272296530693728854503386067454,    -0.7674079606441468267385037317005977222,    0.08645472731882842109826654183371845927,
  -0.7674079606441468267385037317005977222,   -0.48909154272296530693728854503386067454,    0.08645472731882842109826654183371845927,
  0.70267558502048008323600722086411689719,   -0.74544056553282126242431625050523532213,   0.030167355153022877171701180925533825614,
 -0.74544056553282126242431625050523532213,    0.70267558502048008323600722086411689719,   0.030167355153022877171701180925533825614,
 -0.95723501948765882081169097035888157506,   -0.74544056553282126242431625050523532213,   0.030167355153022877171701180925533825614,
 -0.74544056553282126242431625050523532213,   -0.95723501948765882081169097035888157506,   0.030167355153022877171701180925533825614,
 -0.95723501948765882081169097035888157506,    0.70267558502048008323600722086411689719,   0.030167355153022877171701180925533825614,
  0.70267558502048008323600722086411689719,   -0.95723501948765882081169097035888157506,   0.030167355153022877171701180925533825614}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 12 and 13 orthotope
    tempNumIPs = 37;
    nIPMap[cubDataKey(2, 12, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 13, orthotope)] = tempNumIPs;
    // nIP 37
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {                                       0,                                          0,   0.29999334435890964324064310267149902525,
 0.98346133261324557127407851401210794962,                                          0,  0.038128625349852737619351565405684462236,
                                        0,   0.98346133261324557127407851401210794962,  0.038128625349852737619351565405684462236,
-0.98346133261324557127407851401210794962,                                          0,  0.038128625349852737619351565405684462236,
                                        0,  -0.98346133261324557127407851401210794962,  0.038128625349852737619351565405684462236,
 0.63986141836710976174616715579576307077,                                          0,   0.18453546896980721238582475319389847126,
                                        0,   0.63986141836710976174616715579576307077,   0.18453546896980721238582475319389847126,
-0.63986141836710976174616715579576307077,                                          0,   0.18453546896980721238582475319389847126,
                                        0,  -0.63986141836710976174616715579576307077,   0.18453546896980721238582475319389847126,
 0.91877848807971140189839438297320428483,   0.91877848807971140189839438297320428483,  0.039507140647452433689267407341729599229,
 0.91877848807971140189839438297320428483,  -0.91877848807971140189839438297320428483,  0.039507140647452433689267407341729599229,
-0.91877848807971140189839438297320428483,   0.91877848807971140189839438297320428483,  0.039507140647452433689267407341729599229,
-0.91877848807971140189839438297320428483,  -0.91877848807971140189839438297320428483,  0.039507140647452433689267407341729599229,
 0.37962850518674385579275417103548349005,   0.37962850518674385579275417103548349005,   0.23139851940068538741918311064046319298,
 0.37962850518674385579275417103548349005,  -0.37962850518674385579275417103548349005,   0.23139851940068538741918311064046319298,
-0.37962850518674385579275417103548349005,   0.37962850518674385579275417103548349005,   0.23139851940068538741918311064046319298,
-0.37962850518674385579275417103548349005,  -0.37962850518674385579275417103548349005,   0.23139851940068538741918311064046319298,
 0.69955421651335111029920603284399340644,   0.69955421651335111029920603284399340644,   0.13724209671300348677866127359579016247,
 0.69955421651335111029920603284399340644,  -0.69955421651335111029920603284399340644,   0.13724209671300348677866127359579016247,
-0.69955421651335111029920603284399340644,   0.69955421651335111029920603284399340644,   0.13724209671300348677866127359579016247,
-0.69955421651335111029920603284399340644,  -0.69955421651335111029920603284399340644,   0.13724209671300348677866127359579016247,
 0.64359737499661812031249815775600267943,   0.97506888390598351993655206803811094573,  0.033519780050381430279388297868080126418,
 0.97506888390598351993655206803811094573,   0.64359737499661812031249815775600267943,  0.033519780050381430279388297868080126418,
 0.64359737499661812031249815775600267943,  -0.97506888390598351993655206803811094573,  0.033519780050381430279388297868080126418,
-0.97506888390598351993655206803811094573,   0.64359737499661812031249815775600267943,  0.033519780050381430279388297868080126418,
-0.64359737499661812031249815775600267943,   0.97506888390598351993655206803811094573,  0.033519780050381430279388297868080126418,
 0.97506888390598351993655206803811094573,  -0.64359737499661812031249815775600267943,  0.033519780050381430279388297868080126418,
-0.64359737499661812031249815775600267943,  -0.97506888390598351993655206803811094573,  0.033519780050381430279388297868080126418,
-0.97506888390598351993655206803811094573,  -0.64359737499661812031249815775600267943,  0.033519780050381430279388297868080126418,
 0.33353988116478310285167199804895857823,   0.86442760926706191157507267741750728791,   0.11357512636435423536938725920919955133,
 0.86442760926706191157507267741750728791,   0.33353988116478310285167199804895857823,   0.11357512636435423536938725920919955133,
 0.33353988116478310285167199804895857823,  -0.86442760926706191157507267741750728791,   0.11357512636435423536938725920919955133,
-0.86442760926706191157507267741750728791,   0.33353988116478310285167199804895857823,   0.11357512636435423536938725920919955133,
-0.33353988116478310285167199804895857823,   0.86442760926706191157507267741750728791,   0.11357512636435423536938725920919955133,
 0.86442760926706191157507267741750728791,  -0.33353988116478310285167199804895857823,   0.11357512636435423536938725920919955133,
-0.33353988116478310285167199804895857823,  -0.86442760926706191157507267741750728791,   0.11357512636435423536938725920919955133,
-0.86442760926706191157507267741750728791,  -0.33353988116478310285167199804895857823,   0.11357512636435423536938725920919955133}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 13 simplex
    tempNumIPs = 37;
    nIPMap[cubDataKey(2, 13, simplex)] = tempNumIPs;
    // nIP 37
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,  -0.33333333333333333333333333333333333333,   0.13592007317366328856354884936176975902,
-0.02184610709492130019862056181959122469,  -0.95630778581015739960275887636081755062,  0.047988803857789461547421598901919307449,
-0.95630778581015739960275887636081755062,  -0.02184610709492130019862056181959122469,  0.047988803857789461547421598901919307449,
-0.02184610709492130019862056181959122469,  -0.02184610709492130019862056181959122469,  0.047988803857789461547421598901919307449,
-0.55725542741633419869037489058984184448,   0.11451085483266839738074978117968368896,   0.11655697023839996280953416702667961903,
 0.11451085483266839738074978117968368896,  -0.55725542741633419869037489058984184448,   0.11655697023839996280953416702667961903,
-0.55725542741633419869037489058984184448,  -0.55725542741633419869037489058984184448,   0.11655697023839996280953416702667961903,
-0.14611717148039918795837492993725157719,  -0.70776565703920162408325014012549684562,   0.11120393506090665741451493202092293506,
-0.70776565703920162408325014012549684562,  -0.14611717148039918795837492993725157719,   0.11120393506090665741451493202092293506,
-0.14611717148039918795837492993725157719,  -0.14611717148039918795837492993725157719,   0.11120393506090665741451493202092293506,
 -0.9569806377823136322614173729318958348,   0.91396127556462726452283474586379166961,  0.012104674207078343683585600064581648939,
 0.91396127556462726452283474586379166961,   -0.9569806377823136322614173729318958348,  0.012104674207078343683585600064581648939,
 -0.9569806377823136322614173729318958348,   -0.9569806377823136322614173729318958348,  0.012104674207078343683585600064581648939,
-0.82420903393560535081381299499933379777,   0.49701423179990439034603719157741930009,  0.048358079623187638274891491146121518585,
 0.49701423179990439034603719157741930009,  -0.82420903393560535081381299499933379777,  0.048358079623187638274891491146121518585,
-0.67280519786429903953222419657808550232,   0.49701423179990439034603719157741930009,  0.048358079623187638274891491146121518585,
 0.49701423179990439034603719157741930009,  -0.67280519786429903953222419657808550232,  0.048358079623187638274891491146121518585,
-0.67280519786429903953222419657808550232,  -0.82420903393560535081381299499933379777,  0.048358079623187638274891491146121518585,
-0.82420903393560535081381299499933379777,  -0.67280519786429903953222419657808550232,  0.048358079623187638274891491146121518585,
-0.77815591439307320917426090955665095644,   0.72941554059088555060509190179138636667,  0.029930802210331334526491714265806876009,
 0.72941554059088555060509190179138636667,  -0.77815591439307320917426090955665095644,  0.029930802210331334526491714265806876009,
-0.95125962619781234143083099223473541023,   0.72941554059088555060509190179138636667,  0.029930802210331334526491714265806876009,
 0.72941554059088555060509190179138636667,  -0.95125962619781234143083099223473541023,  0.029930802210331334526491714265806876009,
-0.95125962619781234143083099223473541023,  -0.77815591439307320917426090955665095644,  0.029930802210331334526491714265806876009,
-0.77815591439307320917426090955665095644,  -0.95125962619781234143083099223473541023,  0.029930802210331334526491714265806876009,
-0.38311647821576445068305629491750938056,   0.24709199110735114163170870637247317166,  0.069282552281696740931973657021836440891,
 0.24709199110735114163170870637247317166,  -0.38311647821576445068305629491750938056,  0.069282552281696740931973657021836440891,
 -0.8639755128915866909486524114549637911,   0.24709199110735114163170870637247317166,  0.069282552281696740931973657021836440891,
 0.24709199110735114163170870637247317166,   -0.8639755128915866909486524114549637911,  0.069282552281696740931973657021836440891,
 -0.8639755128915866909486524114549637911,  -0.38311647821576445068305629491750938056,  0.069282552281696740931973657021836440891,
-0.38311647821576445068305629491750938056,   -0.8639755128915866909486524114549637911,  0.069282552281696740931973657021836440891,
-0.98974722179523526288134119331774178805,   0.44471558624837593052124026460956809725,  0.019181362007086525445190180332221782772,
 0.44471558624837593052124026460956809725,  -0.98974722179523526288134119331774178805,  0.019181362007086525445190180332221782772,
 -0.4549683644531406676398990712918263092,   0.44471558624837593052124026460956809725,  0.019181362007086525445190180332221782772,
 0.44471558624837593052124026460956809725,   -0.4549683644531406676398990712918263092,  0.019181362007086525445190180332221782772,
 -0.4549683644531406676398990712918263092,  -0.98974722179523526288134119331774178805,  0.019181362007086525445190180332221782772,
-0.98974722179523526288134119331774178805,   -0.4549683644531406676398990712918263092,  0.019181362007086525445190180332221782772}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 14 simplex
    tempNumIPs = 42;
    nIPMap[cubDataKey(2, 14, simplex)] = tempNumIPs;
    // nIP 42
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.64558893517491312608677861906988184783,    0.29117787034982625217355723813976369567,   0.084325177473986035076460874648372259691,
  0.29117787034982625217355723813976369567,   -0.64558893517491312608677861906988184783,   0.084325177473986035076460874648372259691,
 -0.64558893517491312608677861906988184783,   -0.64558893517491312608677861906988184783,   0.084325177473986035076460874648372259691,
 -0.16471056131909215498111835562871312178,   -0.67057887736181569003776328874257375644,   0.065576707088250701282621957477250686472,
 -0.67057887736181569003776328874257375644,   -0.16471056131909215498111835562871312178,   0.065576707088250701282621957477250686472,
 -0.16471056131909215498111835562871312178,   -0.16471056131909215498111835562871312178,   0.065576707088250701282621957477250686472,
 -0.87640023381825479746504234312612842323,    0.75280046763650959493008468625225684647,   0.028867399339553335203419842961306633228,
  0.75280046763650959493008468625225684647,   -0.87640023381825479746504234312612842323,   0.028867399339553335203419842961306633228,
 -0.87640023381825479746504234312612842323,   -0.87640023381825479746504234312612842323,   0.028867399339553335203419842961306633228,
-0.022072179275642722645247959095219528414,   -0.95585564144871455470950408180956094317,   0.043767162738857781281689891926651943066,
 -0.95585564144871455470950408180956094317,  -0.022072179275642722645247959095219528414,   0.043767162738857781281689891926651943066,
-0.022072179275642722645247959095219528414,  -0.022072179275642722645247959095219528414,   0.043767162738857781281689891926651943066,
 -0.45304494338232268049011143347460288761,  -0.093910113235354639019777133050794224782,    0.10354820901458317262956982033279280483,
-0.093910113235354639019777133050794224782,   -0.45304494338232268049011143347460288761,    0.10354820901458317262956982033279280483,
 -0.45304494338232268049011143347460288761,   -0.45304494338232268049011143347460288761,    0.10354820901458317262956982033279280483,
 -0.96121807750259790364349980989094097785,     0.9224361550051958072869996197818819557,   0.009846807204800163363652047018084309531,
   0.9224361550051958072869996197818819557,   -0.96121807750259790364349980989094097785,   0.009846807204800163363652047018084309531,
 -0.96121807750259790364349980989094097785,   -0.96121807750259790364349980989094097785,   0.009846807204800163363652047018084309531,
 -0.40325423572748449405833696388077454282,    0.37396033561617567471725430804062611952,   0.028872616227067680992177383998031606893,
  0.37396033561617567471725430804062611952,   -0.40325423572748449405833696388077454282,   0.028872616227067680992177383998031606893,
 -0.97070609988869118065891734415985157669,    0.37396033561617567471725430804062611952,   0.028872616227067680992177383998031606893,
  0.37396033561617567471725430804062611952,   -0.97070609988869118065891734415985157669,   0.028872616227067680992177383998031606893,
 -0.97070609988869118065891734415985157669,   -0.40325423572748449405833696388077454282,   0.028872616227067680992177383998031606893,
 -0.40325423572748449405833696388077454282,   -0.97070609988869118065891734415985157669,   0.028872616227067680992177383998031606893,
 -0.88575048519270412192864575156217057937,    0.54121710954999296517806654833485591265,   0.049331506425127347925750490367272456636,
  0.54121710954999296517806654833485591265,   -0.88575048519270412192864575156217057937,   0.049331506425127347925750490367272456636,
 -0.65546662435728884324942079677268533329,    0.54121710954999296517806654833485591265,   0.049331506425127347925750490367272456636,
  0.54121710954999296517806654833485591265,   -0.65546662435728884324942079677268533329,   0.049331506425127347925750490367272456636,
 -0.65546662435728884324942079677268533329,   -0.88575048519270412192864575156217057937,   0.049331506425127347925750490367272456636,
 -0.88575048519270412192864575156217057937,   -0.65546662435728884324942079677268533329,   0.049331506425127347925750490367272456636,
 -0.32627708040730999651188960582214922188,    0.14044458169336634699539242672470851928,   0.077143021574121366456978055620821718469,
  0.14044458169336634699539242672470851928,   -0.32627708040730999651188960582214922188,   0.077143021574121366456978055620821718469,
  -0.8141675012860563504835028209025592974,    0.14044458169336634699539242672470851928,   0.077143021574121366456978055620821718469,
  0.14044458169336634699539242672470851928,    -0.8141675012860563504835028209025592974,   0.077143021574121366456978055620821718469,
  -0.8141675012860563504835028209025592974,   -0.32627708040730999651188960582214922188,   0.077143021574121366456978055620821718469,
 -0.32627708040730999651188960582214922188,    -0.8141675012860563504835028209025592974,   0.077143021574121366456978055620821718469,
 -0.99746333813425594982550719780901461465,    0.75951434274034225902914327394920365201,   0.010020457677001343539720186164978232925,
  0.75951434274034225902914327394920365201,   -0.99746333813425594982550719780901461465,   0.010020457677001343539720186164978232925,
 -0.76205100460608630920363607614018903735,    0.75951434274034225902914327394920365201,   0.010020457677001343539720186164978232925,
  0.75951434274034225902914327394920365201,   -0.76205100460608630920363607614018903735,   0.010020457677001343539720186164978232925,
 -0.76205100460608630920363607614018903735,   -0.99746333813425594982550719780901461465,   0.010020457677001343539720186164978232925,
 -0.99746333813425594982550719780901461465,   -0.76205100460608630920363607614018903735,   0.010020457677001343539720186164978232925}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 14 and 15 orthotope
    tempNumIPs = 48;
    nIPMap[cubDataKey(2, 14, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 15, orthotope)] = tempNumIPs;
    // nIP 48
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0.79809898789700709024555173503596148462,                                           0,    0.11493292726352536781596028932541766228,
                                         0,    0.79809898789700709024555173503596148462,    0.11493292726352536781596028932541766228,
 -0.79809898789700709024555173503596148462,                                           0,    0.11493292726352536781596028932541766228,
                                         0,   -0.79809898789700709024555173503596148462,    0.11493292726352536781596028932541766228,
  0.30385302639845975610268178690962703058,                                           0,    0.18168575896727184957040458577073367484,
                                         0,    0.30385302639845975610268178690962703058,    0.18168575896727184957040458577073367484,
 -0.30385302639845975610268178690962703058,                                           0,    0.18168575896727184957040458577073367484,
                                         0,   -0.30385302639845975610268178690962703058,    0.18168575896727184957040458577073367484,
  0.88242227010688535992041188667767184542,    0.88242227010688535992041188667767184542,   0.041238483788765811839027639617216599853,
  0.88242227010688535992041188667767184542,   -0.88242227010688535992041188667767184542,   0.041238483788765811839027639617216599853,
 -0.88242227010688535992041188667767184542,    0.88242227010688535992041188667767184542,   0.041238483788765811839027639617216599853,
 -0.88242227010688535992041188667767184542,   -0.88242227010688535992041188667767184542,   0.041238483788765811839027639617216599853,
  0.97778979953990269742424689352375866785,    0.97778979953990269742424689352375866785,  0.0059337464858392213147614197893498681047,
  0.97778979953990269742424689352375866785,   -0.97778979953990269742424689352375866785,  0.0059337464858392213147614197893498681047,
 -0.97778979953990269742424689352375866785,    0.97778979953990269742424689352375866785,  0.0059337464858392213147614197893498681047,
 -0.97778979953990269742424689352375866785,   -0.97778979953990269742424689352375866785,  0.0059337464858392213147614197893498681047,
  0.80871213581950356724002823901416977471,    0.56721357339659070317865169083903012406,    0.10114914347743242280761079470300175928,
  0.56721357339659070317865169083903012406,    0.80871213581950356724002823901416977471,    0.10114914347743242280761079470300175928,
  0.80871213581950356724002823901416977471,   -0.56721357339659070317865169083903012406,    0.10114914347743242280761079470300175928,
 -0.56721357339659070317865169083903012406,    0.80871213581950356724002823901416977471,    0.10114914347743242280761079470300175928,
 -0.80871213581950356724002823901416977471,    0.56721357339659070317865169083903012406,    0.10114914347743242280761079470300175928,
  0.56721357339659070317865169083903012406,   -0.80871213581950356724002823901416977471,    0.10114914347743242280761079470300175928,
 -0.80871213581950356724002823901416977471,   -0.56721357339659070317865169083903012406,    0.10114914347743242280761079470300175928,
 -0.56721357339659070317865169083903012406,   -0.80871213581950356724002823901416977471,    0.10114914347743242280761079470300175928,
  0.30465479903705254117098598695333881061,    0.57873619403580651026013939264701051891,    0.14783672163288121847905618360636339726,
  0.57873619403580651026013939264701051891,    0.30465479903705254117098598695333881061,    0.14783672163288121847905618360636339726,
  0.30465479903705254117098598695333881061,   -0.57873619403580651026013939264701051891,    0.14783672163288121847905618360636339726,
 -0.57873619403580651026013939264701051891,    0.30465479903705254117098598695333881061,    0.14783672163288121847905618360636339726,
 -0.30465479903705254117098598695333881061,    0.57873619403580651026013939264701051891,    0.14783672163288121847905618360636339726,
  0.57873619403580651026013939264701051891,   -0.30465479903705254117098598695333881061,    0.14783672163288121847905618360636339726,
 -0.30465479903705254117098598695333881061,   -0.57873619403580651026013939264701051891,    0.14783672163288121847905618360636339726,
 -0.57873619403580651026013939264701051891,   -0.30465479903705254117098598695333881061,    0.14783672163288121847905618360636339726,
   0.9805048437245319615400269877302555585,     0.6974636319909670187938904551822360596,   0.023273194144673214602574891933463059923,
   0.6974636319909670187938904551822360596,     0.9805048437245319615400269877302555585,   0.023273194144673214602574891933463059923,
   0.9805048437245319615400269877302555585,    -0.6974636319909670187938904551822360596,   0.023273194144673214602574891933463059923,
  -0.6974636319909670187938904551822360596,     0.9805048437245319615400269877302555585,   0.023273194144673214602574891933463059923,
  -0.9805048437245319615400269877302555585,     0.6974636319909670187938904551822360596,   0.023273194144673214602574891933463059923,
   0.6974636319909670187938904551822360596,    -0.9805048437245319615400269877302555585,   0.023273194144673214602574891933463059923,
  -0.9805048437245319615400269877302555585,    -0.6974636319909670187938904551822360596,   0.023273194144673214602574891933463059923,
  -0.6974636319909670187938904551822360596,    -0.9805048437245319615400269877302555585,   0.023273194144673214602574891933463059923,
  0.26484415587231617710014274826830856062,    0.95574251980951170852387507604700262081,   0.055845482492312018840681162505812880996,
  0.95574251980951170852387507604700262081,    0.26484415587231617710014274826830856062,   0.055845482492312018840681162505812880996,
  0.26484415587231617710014274826830856062,   -0.95574251980951170852387507604700262081,   0.055845482492312018840681162505812880996,
 -0.95574251980951170852387507604700262081,    0.26484415587231617710014274826830856062,   0.055845482492312018840681162505812880996,
 -0.26484415587231617710014274826830856062,    0.95574251980951170852387507604700262081,   0.055845482492312018840681162505812880996,
  0.95574251980951170852387507604700262081,   -0.26484415587231617710014274826830856062,   0.055845482492312018840681162505812880996,
 -0.26484415587231617710014274826830856062,   -0.95574251980951170852387507604700262081,   0.055845482492312018840681162505812880996,
 -0.95574251980951170852387507604700262081,   -0.26484415587231617710014274826830856062,   0.055845482492312018840681162505812880996}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 15 simplex
    tempNumIPs = 49;
    nIPMap[cubDataKey(2, 15, simplex)] = tempNumIPs;
    // nIP 49
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,    0.08867077476436815077519868131406480795,
 -0.18927557173204904426648777856049828207,   -0.62144885653590191146702444287900343586,   0.085427563142921134082725838745668360696,
 -0.62144885653590191146702444287900343586,   -0.18927557173204904426648777856049828207,   0.085427563142921134082725838745668360696,
 -0.18927557173204904426648777856049828207,   -0.18927557173204904426648777856049828207,   0.085427563142921134082725838745668360696,
 -0.85965289420002796503128049571314279073,    0.71930578840005593006256099142628558146,    0.03288947512525032753736001737701825892,
  0.71930578840005593006256099142628558146,   -0.85965289420002796503128049571314279073,    0.03288947512525032753736001737701825892,
 -0.85965289420002796503128049571314279073,   -0.85965289420002796503128049571314279073,    0.03288947512525032753736001737701825892,
-0.051658637123960416042618609451788188253,   -0.89668272575207916791476278109642362349,   0.034792296001526829115453101768509506658,
 -0.89668272575207916791476278109642362349,  -0.051658637123960416042618609451788188253,   0.034792296001526829115453101768509506658,
-0.051658637123960416042618609451788188253,  -0.051658637123960416042618609451788188253,   0.034792296001526829115453101768509506658,
 -0.54724257315930080341854883783076139455,   0.094485146318601606837097675661522789095,   0.093566723457419264712537789561548478374,
 0.094485146318601606837097675661522789095,   -0.54724257315930080341854883783076139455,   0.093566723457419264712537789561548478374,
 -0.54724257315930080341854883783076139455,   -0.54724257315930080341854883783076139455,   0.093566723457419264712537789561548478374,
-0.010006086461747631524400704632178282105,   -0.97998782707650473695119859073564343579,   0.019147692364920172413138507253187117414,
 -0.97998782707650473695119859073564343579,  -0.010006086461747631524400704632178282105,   0.019147692364920172413138507253187117414,
-0.010006086461747631524400704632178282105,  -0.010006086461747631524400704632178282105,   0.019147692364920172413138507253187117414,
 -0.96837654749802270888209653243226493969,    0.93675309499604541776419306486452987938,  0.0059215492758107511385526052767298787739,
  0.93675309499604541776419306486452987938,   -0.96837654749802270888209653243226493969,  0.0059215492758107511385526052767298787739,
 -0.96837654749802270888209653243226493969,   -0.96837654749802270888209653243226493969,  0.0059215492758107511385526052767298787739,
  -0.9632477752286378160779112088372081634,    0.33395128960373617989264707840742420014,   0.031205145661151929503549831657000846318,
  0.33395128960373617989264707840742420014,    -0.9632477752286378160779112088372081634,   0.031205145661151929503549831657000846318,
 -0.37070351437509836381473586957021603674,    0.33395128960373617989264707840742420014,   0.031205145661151929503549831657000846318,
  0.33395128960373617989264707840742420014,   -0.37070351437509836381473586957021603674,   0.031205145661151929503549831657000846318,
 -0.37070351437509836381473586957021603674,    -0.9632477752286378160779112088372081634,   0.031205145661151929503549831657000846318,
  -0.9632477752286378160779112088372081634,   -0.37070351437509836381473586957021603674,   0.031205145661151929503549831657000846318,
 -0.98172152592538320992689636501802696381,    0.83982431545247212980572119059358474792,  0.0080597067440361992367286698217130156564,
  0.83982431545247212980572119059358474792,   -0.98172152592538320992689636501802696381,  0.0080597067440361992367286698217130156564,
 -0.85810278952708891987882482557555778411,    0.83982431545247212980572119059358474792,  0.0080597067440361992367286698217130156564,
  0.83982431545247212980572119059358474792,   -0.85810278952708891987882482557555778411,  0.0080597067440361992367286698217130156564,
 -0.85810278952708891987882482557555778411,   -0.98172152592538320992689636501802696381,  0.0080597067440361992367286698217130156564,
 -0.98172152592538320992689636501802696381,   -0.85810278952708891987882482557555778411,  0.0080597067440361992367286698217130156564,
  -0.6189288210472121331386195006280717197,    0.43044471386290140376795873756746088333,   0.057441173850402680461331283343301712008,
  0.43044471386290140376795873756746088333,    -0.6189288210472121331386195006280717197,   0.057441173850402680461331283343301712008,
 -0.81151589281568927062933923693938916363,    0.43044471386290140376795873756746088333,   0.057441173850402680461331283343301712008,
  0.43044471386290140376795873756746088333,   -0.81151589281568927062933923693938916363,   0.057441173850402680461331283343301712008,
 -0.81151589281568927062933923693938916363,    -0.6189288210472121331386195006280717197,   0.057441173850402680461331283343301712008,
  -0.6189288210472121331386195006280717197,   -0.81151589281568927062933923693938916363,   0.057441173850402680461331283343301712008,
 -0.66386270955517128907044482283274107879,    0.62658528209883853982451662602671593723,   0.023345242363151690679154011088673171842,
  0.62658528209883853982451662602671593723,   -0.66386270955517128907044482283274107879,   0.023345242363151690679154011088673171842,
 -0.96272257254366725075407180319397485844,    0.62658528209883853982451662602671593723,   0.023345242363151690679154011088673171842,
  0.62658528209883853982451662602671593723,   -0.96272257254366725075407180319397485844,   0.023345242363151690679154011088673171842,
 -0.96272257254366725075407180319397485844,   -0.66386270955517128907044482283274107879,   0.023345242363151690679154011088673171842,
 -0.66386270955517128907044482283274107879,   -0.96272257254366725075407180319397485844,   0.023345242363151690679154011088673171842,
 -0.32209877704944562157081399065909603587,    0.13050532975422842890370386873351257506,   0.062630952569938568823485827212302985766,
  0.13050532975422842890370386873351257506,   -0.32209877704944562157081399065909603587,   0.062630952569938568823485827212302985766,
 -0.80840655270478280733288987807441653919,    0.13050532975422842890370386873351257506,   0.062630952569938568823485827212302985766,
  0.13050532975422842890370386873351257506,   -0.80840655270478280733288987807441653919,   0.062630952569938568823485827212302985766,
 -0.80840655270478280733288987807441653919,   -0.32209877704944562157081399065909603587,   0.062630952569938568823485827212302985766,
 -0.32209877704944562157081399065909603587,   -0.80840655270478280733288987807441653919,   0.062630952569938568823485827212302985766}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 16 simplex
    tempNumIPs = 55;
    nIPMap[cubDataKey(2, 16, simplex)] = tempNumIPs;
    // nIP 55
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,      -0.33333333333333333333333333333333333333,      0.090529132147637585798682804751497905871,
      -0.508019859065716531056564440465839641,      0.016039718131433062113128880931679281991,      0.082185846287397899082859837743846349856,
    0.016039718131433062113128880931679281991,        -0.508019859065716531056564440465839641,      0.082185846287397899082859837743846349856,
      -0.508019859065716531056564440465839641,        -0.508019859065716531056564440465839641,      0.082185846287397899082859837743846349856,
     -0.1688302062291589735580386476319745286,       -0.6623395875416820528839227047360509428,      0.081423666624850718504727604609892231433,
     -0.6623395875416820528839227047360509428,       -0.1688302062291589735580386476319745286,      0.081423666624850718504727604609892231433,
     -0.1688302062291589735580386476319745286,       -0.1688302062291589735580386476319745286,      0.081423666624850718504727604609892231433,
    -0.82928886682659935648935022895726939779,       0.65857773365319871297870045791453879558,      0.029563269380448807960364375206438726749,
     0.65857773365319871297870045791453879558,      -0.82928886682659935648935022895726939779,      0.029563269380448807960364375206438726749,
    -0.82928886682659935648935022895726939779,      -0.82928886682659935648935022895726939779,      0.029563269380448807960364375206438726749,
    -0.67616271161745754624623346597531512883,       0.35232542323491509249246693195063025766,      0.058836819397976195721222704445565411515,
     0.35232542323491509249246693195063025766,      -0.67616271161745754624623346597531512883,      0.058836819397976195721222704445565411515,
    -0.67616271161745754624623346597531512883,      -0.67616271161745754624623346597531512883,      0.058836819397976195721222704445565411515,
-5.3815663675432935238750288859902834282e-109,                                             -1,     0.0088370926243011378507684635783092181714,
                                           -1,  -5.3815663675432935238750288859902834282e-109,     0.0088370926243011378507684635783092181714,
-5.3815663675432935238750288859902834282e-109,  -5.3815663675432935238750288859902834282e-109,     0.0088370926243011378507684635783092181714,
   -0.049438544908115791823484418930976089946,      -0.90112291018376841635303116213804782011,      0.051948666596554316640604583694066880977,
    -0.90112291018376841635303116213804782011,     -0.049438544908115791823484418930976089946,      0.051948666596554316640604583694066880977,
   -0.049438544908115791823484418930976089946,     -0.049438544908115791823484418930976089946,      0.051948666596554316640604583694066880977,
    -0.89048965017059377421655257720683258971,       0.50834012288953561570997807325829055003,      0.037876544928831397939092427613418818028,
     0.50834012288953561570997807325829055003,      -0.89048965017059377421655257720683258971,      0.037876544928831397939092427613418818028,
    -0.61785047271894184149342549605145796032,       0.50834012288953561570997807325829055003,      0.037876544928831397939092427613418818028,
     0.50834012288953561570997807325829055003,      -0.61785047271894184149342549605145796032,      0.037876544928831397939092427613418818028,
    -0.61785047271894184149342549605145796032,      -0.89048965017059377421655257720683258971,      0.037876544928831397939092427613418818028,
    -0.89048965017059377421655257720683258971,      -0.61785047271894184149342549605145796032,      0.037876544928831397939092427613418818028,
    -0.95359314446237258192502291745679921363,       0.93648873606191735836919518841625197884,     0.0033089334296700963847916278962044308229,
     0.93648873606191735836919518841625197884,      -0.95359314446237258192502291745679921363,     0.0033089334296700963847916278962044308229,
    -0.98289559159954477644417227095945276521,       0.93648873606191735836919518841625197884,     0.0033089334296700963847916278962044308229,
     0.93648873606191735836919518841625197884,      -0.98289559159954477644417227095945276521,     0.0033089334296700963847916278962044308229,
    -0.98289559159954477644417227095945276521,      -0.95359314446237258192502291745679921363,     0.0033089334296700963847916278962044308229,
    -0.95359314446237258192502291745679921363,      -0.98289559159954477644417227095945276521,     0.0033089334296700963847916278962044308229,
    -0.96213644343918817003656945894522367429,       0.29860739649089281625925299156381596948,      0.030017203568571612165136964521714949941,
     0.29860739649089281625925299156381596948,      -0.96213644343918817003656945894522367429,      0.030017203568571612165136964521714949941,
     -0.3364709530517046462226835326185922952,       0.29860739649089281625925299156381596948,      0.030017203568571612165136964521714949941,
     0.29860739649089281625925299156381596948,       -0.3364709530517046462226835326185922952,      0.030017203568571612165136964521714949941,
     -0.3364709530517046462226835326185922952,      -0.96213644343918817003656945894522367429,      0.030017203568571612165136964521714949941,
    -0.96213644343918817003656945894522367429,       -0.3364709530517046462226835326185922952,      0.030017203568571612165136964521714949941,
    -0.96193974051260510241812957568503652064,       0.80054740654085925227401383438774365525,      0.015895187866784997621199196218108093516,
     0.80054740654085925227401383438774365525,      -0.96193974051260510241812957568503652064,      0.015895187866784997621199196218108093516,
    -0.83860766602825414985588425870270713461,       0.80054740654085925227401383438774365525,      0.015895187866784997621199196218108093516,
     0.80054740654085925227401383438774365525,      -0.83860766602825414985588425870270713461,      0.015895187866784997621199196218108093516,
    -0.83860766602825414985588425870270713461,      -0.96193974051260510241812957568503652064,      0.015895187866784997621199196218108093516,
    -0.96193974051260510241812957568503652064,      -0.83860766602825414985588425870270713461,      0.015895187866784997621199196218108093516,
    -0.79478761952120380739655029677085645407,       0.17829768112849583758057340125685564327,      0.063967220158740137569628123287866229484,
     0.17829768112849583758057340125685564327,      -0.79478761952120380739655029677085645407,      0.063967220158740137569628123287866229484,
     -0.3835100616072920301840231044859991892,       0.17829768112849583758057340125685564327,      0.063967220158740137569628123287866229484,
     0.17829768112849583758057340125685564327,       -0.3835100616072920301840231044859991892,      0.063967220158740137569628123287866229484,
     -0.3835100616072920301840231044859991892,      -0.79478761952120380739655029677085645407,      0.063967220158740137569628123287866229484,
    -0.79478761952120380739655029677085645407,       -0.3835100616072920301840231044859991892,      0.063967220158740137569628123287866229484,
      -0.988127299966355460185036064778683968,       0.61324373499879128490286065261412937757,      0.010782374233697622806764075031711751212,
     0.61324373499879128490286065261412937757,        -0.988127299966355460185036064778683968,      0.010782374233697622806764075031711751212,
    -0.62511643503243582471782458783544540958,       0.61324373499879128490286065261412937757,      0.010782374233697622806764075031711751212,
     0.61324373499879128490286065261412937757,      -0.62511643503243582471782458783544540958,      0.010782374233697622806764075031711751212,
    -0.62511643503243582471782458783544540958,        -0.988127299966355460185036064778683968,      0.010782374233697622806764075031711751212,
      -0.988127299966355460185036064778683968,      -0.62511643503243582471782458783544540958,      0.010782374233697622806764075031711751212}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 16 and 17 orthotope
    tempNumIPs = 60;
    nIPMap[cubDataKey(2, 16, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 17, orthotope)] = tempNumIPs;
    // nIP 60
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0.54504298720914838787010174782667747048,                                           0,    0.13880437768721992480048871346267772478,
                                         0,    0.54504298720914838787010174782667747048,    0.13880437768721992480048871346267772478,
 -0.54504298720914838787010174782667747048,                                           0,    0.13880437768721992480048871346267772478,
                                         0,   -0.54504298720914838787010174782667747048,    0.13880437768721992480048871346267772478,
  0.91749718821915320817120277758673277115,                                           0,   0.065347823807464912859066688973533173055,
                                         0,    0.91749718821915320817120277758673277115,   0.065347823807464912859066688973533173055,
 -0.91749718821915320817120277758673277115,                                           0,   0.065347823807464912859066688973533173055,
                                         0,   -0.91749718821915320817120277758673277115,   0.065347823807464912859066688973533173055,
  0.19255424261404688054827243228978233862,    0.19255424261404688054827243228978233862,    0.14285641962219959571166270516073482911,
  0.19255424261404688054827243228978233862,   -0.19255424261404688054827243228978233862,    0.14285641962219959571166270516073482911,
 -0.19255424261404688054827243228978233862,    0.19255424261404688054827243228978233862,    0.14285641962219959571166270516073482911,
 -0.19255424261404688054827243228978233862,   -0.19255424261404688054827243228978233862,    0.14285641962219959571166270516073482911,
  0.87132868469544327271341598059215273696,    0.87132868469544327271341598059215273696,   0.036715908485251512466108823252963750875,
  0.87132868469544327271341598059215273696,   -0.87132868469544327271341598059215273696,   0.036715908485251512466108823252963750875,
 -0.87132868469544327271341598059215273696,    0.87132868469544327271341598059215273696,   0.036715908485251512466108823252963750875,
 -0.87132868469544327271341598059215273696,   -0.87132868469544327271341598059215273696,   0.036715908485251512466108823252963750875,
   0.6943880900895507516564691872484594143,     0.6943880900895507516564691872484594143,   0.086255699638143804540935141271187587589,
   0.6943880900895507516564691872484594143,    -0.6943880900895507516564691872484594143,   0.086255699638143804540935141271187587589,
  -0.6943880900895507516564691872484594143,     0.6943880900895507516564691872484594143,   0.086255699638143804540935141271187587589,
  -0.6943880900895507516564691872484594143,    -0.6943880900895507516564691872484594143,   0.086255699638143804540935141271187587589,
  0.45768299288315551760664997201644894176,    0.45768299288315551760664997201644894176,    0.13455627321191692494394965878086112148,
  0.45768299288315551760664997201644894176,   -0.45768299288315551760664997201644894176,    0.13455627321191692494394965878086112148,
 -0.45768299288315551760664997201644894176,    0.45768299288315551760664997201644894176,    0.13455627321191692494394965878086112148,
 -0.45768299288315551760664997201644894176,   -0.45768299288315551760664997201644894176,    0.13455627321191692494394965878086112148,
  0.96920425609150873396809064886215595774,    0.96920425609150873396809064886215595774,  0.0076556334149099368494439596667976876552,
  0.96920425609150873396809064886215595774,   -0.96920425609150873396809064886215595774,  0.0076556334149099368494439596667976876552,
 -0.96920425609150873396809064886215595774,    0.96920425609150873396809064886215595774,  0.0076556334149099368494439596667976876552,
 -0.96920425609150873396809064886215595774,   -0.96920425609150873396809064886215595774,  0.0076556334149099368494439596667976876552,
  0.76208827608019827579179765943659247047,    0.27947776671519215001028373596083402351,    0.10317672478493534779017919278157584697,
  0.27947776671519215001028373596083402351,    0.76208827608019827579179765943659247047,    0.10317672478493534779017919278157584697,
  0.76208827608019827579179765943659247047,   -0.27947776671519215001028373596083402351,    0.10317672478493534779017919278157584697,
 -0.27947776671519215001028373596083402351,    0.76208827608019827579179765943659247047,    0.10317672478493534779017919278157584697,
 -0.76208827608019827579179765943659247047,    0.27947776671519215001028373596083402351,    0.10317672478493534779017919278157584697,
  0.27947776671519215001028373596083402351,   -0.76208827608019827579179765943659247047,    0.10317672478493534779017919278157584697,
 -0.76208827608019827579179765943659247047,   -0.27947776671519215001028373596083402351,    0.10317672478493534779017919278157584697,
 -0.27947776671519215001028373596083402351,   -0.76208827608019827579179765943659247047,    0.10317672478493534779017919278157584697,
   0.9073943521982236707040505947884667131,    0.54689867041431010640112807851437665884,   0.054826400237291440916675995759834604896,
  0.54689867041431010640112807851437665884,     0.9073943521982236707040505947884667131,   0.054826400237291440916675995759834604896,
   0.9073943521982236707040505947884667131,   -0.54689867041431010640112807851437665884,   0.054826400237291440916675995759834604896,
 -0.54689867041431010640112807851437665884,     0.9073943521982236707040505947884667131,   0.054826400237291440916675995759834604896,
  -0.9073943521982236707040505947884667131,    0.54689867041431010640112807851437665884,   0.054826400237291440916675995759834604896,
  0.54689867041431010640112807851437665884,    -0.9073943521982236707040505947884667131,   0.054826400237291440916675995759834604896,
  -0.9073943521982236707040505947884667131,   -0.54689867041431010640112807851437665884,   0.054826400237291440916675995759834604896,
 -0.54689867041431010640112807851437665884,    -0.9073943521982236707040505947884667131,   0.054826400237291440916675995759834604896,
  0.76124866643908275965614948240512547767,    0.98378797150154955261318330664974572171,   0.015119810388256861898635715383645346634,
  0.98378797150154955261318330664974572171,    0.76124866643908275965614948240512547767,   0.015119810388256861898635715383645346634,
  0.76124866643908275965614948240512547767,   -0.98378797150154955261318330664974572171,   0.015119810388256861898635715383645346634,
 -0.98378797150154955261318330664974572171,    0.76124866643908275965614948240512547767,   0.015119810388256861898635715383645346634,
 -0.76124866643908275965614948240512547767,    0.98378797150154955261318330664974572171,   0.015119810388256861898635715383645346634,
  0.98378797150154955261318330664974572171,   -0.76124866643908275965614948240512547767,   0.015119810388256861898635715383645346634,
 -0.76124866643908275965614948240512547767,   -0.98378797150154955261318330664974572171,   0.015119810388256861898635715383645346634,
 -0.98378797150154955261318330664974572171,   -0.76124866643908275965614948240512547767,   0.015119810388256861898635715383645346634,
   0.9870947070447678254568263327892903071,    0.30280748178886139816768835788736238796,   0.020780996655963043308681250790566264232,
  0.30280748178886139816768835788736238796,     0.9870947070447678254568263327892903071,   0.020780996655963043308681250790566264232,
   0.9870947070447678254568263327892903071,   -0.30280748178886139816768835788736238796,   0.020780996655963043308681250790566264232,
 -0.30280748178886139816768835788736238796,     0.9870947070447678254568263327892903071,   0.020780996655963043308681250790566264232,
  -0.9870947070447678254568263327892903071,    0.30280748178886139816768835788736238796,   0.020780996655963043308681250790566264232,
  0.30280748178886139816768835788736238796,    -0.9870947070447678254568263327892903071,   0.020780996655963043308681250790566264232,
  -0.9870947070447678254568263327892903071,   -0.30280748178886139816768835788736238796,   0.020780996655963043308681250790566264232,
 -0.30280748178886139816768835788736238796,    -0.9870947070447678254568263327892903071,   0.020780996655963043308681250790566264232}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 17 simplex
    tempNumIPs = 60;
    nIPMap[cubDataKey(2, 17, simplex)] = tempNumIPs;
    // nIP 60
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.16579311127680159678975416725600076679,   -0.66841377744639680642049166548799846641,   0.054621853056204215047557695108059484685,
 -0.66841377744639680642049166548799846641,   -0.16579311127680159678975416725600076679,   0.054621853056204215047557695108059484685,
 -0.16579311127680159678975416725600076679,   -0.16579311127680159678975416725600076679,   0.054621853056204215047557695108059484685,
 -0.97048901667849209319369164139840939816,    0.94097803335698418638738328279681879631,  0.0055477751552752843092249851846592899136,
  0.94097803335698418638738328279681879631,   -0.97048901667849209319369164139840939816,  0.0055477751552752843092249851846592899136,
 -0.97048901667849209319369164139840939816,   -0.97048901667849209319369164139840939816,  0.0055477751552752843092249851846592899136,
-0.068804256762219396208216818306124542516,   -0.86239148647556120758356636338775091497,   0.050038901900994715594115060352765454535,
 -0.86239148647556120758356636338775091497,  -0.068804256762219396208216818306124542516,   0.050038901900994715594115060352765454535,
-0.068804256762219396208216818306124542516,  -0.068804256762219396208216818306124542516,   0.050038901900994715594115060352765454535,
 -0.63928376746725875962831385554123912213,    0.27856753493451751925662771108247824427,   0.052625261176035969912068830347159305457,
  0.27856753493451751925662771108247824427,   -0.63928376746725875962831385554123912213,   0.052625261176035969912068830347159305457,
 -0.63928376746725875962831385554123912213,   -0.63928376746725875962831385554123912213,   0.052625261176035969912068830347159305457,
 -0.86669187304080614048646650320300230601,    0.73338374608161228097293300640600461202,   0.024918001604610884190500085730827367033,
  0.73338374608161228097293300640600461202,   -0.86669187304080614048646650320300230601,   0.024918001604610884190500085730827367033,
 -0.86669187304080614048646650320300230601,   -0.86669187304080614048646650320300230601,   0.024918001604610884190500085730827367033,
 -0.42858699512682674399632424251061595821,   -0.14282600974634651200735151497876808358,   0.075432474305590560032857611935910459517,
 -0.14282600974634651200735151497876808358,   -0.42858699512682674399632424251061595821,   0.075432474305590560032857611935910459517,
 -0.42858699512682674399632424251061595821,   -0.42858699512682674399632424251061595821,   0.075432474305590560032857611935910459517,
 -0.96796471527576140631518537573347042693,     0.6495801403301760488685812163616402457,   0.015956600411859186648565244290131049809,
   0.6495801403301760488685812163616402457,   -0.96796471527576140631518537573347042693,   0.015956600411859186648565244290131049809,
 -0.68161542505441464255339584062816981878,     0.6495801403301760488685812163616402457,   0.015956600411859186648565244290131049809,
   0.6495801403301760488685812163616402457,   -0.68161542505441464255339584062816981878,   0.015956600411859186648565244290131049809,
 -0.68161542505441464255339584062816981878,   -0.96796471527576140631518537573347042693,   0.015956600411859186648565244290131049809,
 -0.96796471527576140631518537573347042693,   -0.68161542505441464255339584062816981878,   0.015956600411859186648565244290131049809,
 -0.38743681650762691692938101229025325485,    0.25273806077290452413744460559562465037,   0.044975545093382132884644311560338005954,
  0.25273806077290452413744460559562465037,   -0.38743681650762691692938101229025325485,   0.044975545093382132884644311560338005954,
 -0.86530124426527760720806359330537139552,    0.25273806077290452413744460559562465037,   0.044975545093382132884644311560338005954,
  0.25273806077290452413744460559562465037,   -0.86530124426527760720806359330537139552,   0.044975545093382132884644311560338005954,
 -0.86530124426527760720806359330537139552,   -0.38743681650762691692938101229025325485,   0.044975545093382132884644311560338005954,
 -0.38743681650762691692938101229025325485,   -0.86530124426527760720806359330537139552,   0.044975545093382132884644311560338005954,
 -0.97354065447982621338596657330178882606,    0.14258973588936810370129822339191719484,   0.020796879911679072927257130756176547979,
  0.14258973588936810370129822339191719484,   -0.97354065447982621338596657330178882606,   0.020796879911679072927257130756176547979,
 -0.16904908140954189031533165009012836878,    0.14258973588936810370129822339191719484,   0.020796879911679072927257130756176547979,
  0.14258973588936810370129822339191719484,   -0.16904908140954189031533165009012836878,   0.020796879911679072927257130756176547979,
 -0.16904908140954189031533165009012836878,   -0.97354065447982621338596657330178882606,   0.020796879911679072927257130756176547979,
 -0.97354065447982621338596657330178882606,   -0.16904908140954189031533165009012836878,   0.020796879911679072927257130756176547979,
 -0.84391531886343515165733825896137376978,    0.50647029187291623779285479822871979869,   0.041115796640909034992011659162843848261,
  0.50647029187291623779285479822871979869,   -0.84391531886343515165733825896137376978,   0.041115796640909034992011659162843848261,
 -0.66255497300948108613551653926734602891,    0.50647029187291623779285479822871979869,   0.041115796640909034992011659162843848261,
  0.50647029187291623779285479822871979869,   -0.66255497300948108613551653926734602891,   0.041115796640909034992011659162843848261,
 -0.66255497300948108613551653926734602891,   -0.84391531886343515165733825896137376978,   0.041115796640909034992011659162843848261,
 -0.84391531886343515165733825896137376978,   -0.66255497300948108613551653926734602891,   0.041115796640909034992011659162843848261,
  -0.9737282583319946100667927211400606547,    0.43014451822128492303184924075272604399,   0.017384429002002383137547850880991473289,
  0.43014451822128492303184924075272604399,    -0.9737282583319946100667927211400606547,   0.017384429002002383137547850880991473289,
 -0.45641625988929031296505651961266538929,    0.43014451822128492303184924075272604399,   0.017384429002002383137547850880991473289,
  0.43014451822128492303184924075272604399,   -0.45641625988929031296505651961266538929,   0.017384429002002383137547850880991473289,
 -0.45641625988929031296505651961266538929,    -0.9737282583319946100667927211400606547,   0.017384429002002383137547850880991473289,
  -0.9737282583319946100667927211400606547,   -0.45641625988929031296505651961266538929,   0.017384429002002383137547850880991473289,
  -0.9768496481936387693039157631967710822,    0.83183870659563391675486447827828667924,  0.0091686968034717336858992147415592016162,
  0.83183870659563391675486447827828667924,    -0.9768496481936387693039157631967710822,  0.0091686968034717336858992147415592016162,
 -0.85498905840199514745094871508151559704,    0.83183870659563391675486447827828667924,  0.0091686968034717336858992147415592016162,
  0.83183870659563391675486447827828667924,   -0.85498905840199514745094871508151559704,  0.0091686968034717336858992147415592016162,
 -0.85498905840199514745094871508151559704,    -0.9768496481936387693039157631967710822,  0.0091686968034717336858992147415592016162,
  -0.9768496481936387693039157631967710822,   -0.85498905840199514745094871508151559704,  0.0091686968034717336858992147415592016162,
 -0.68498904414626018994477816895259857135,   0.086551159192319548209296799380484667575,   0.052343251870673974514245787611602525856,
 0.086551159192319548209296799380484667575,   -0.68498904414626018994477816895259857135,   0.052343251870673974514245787611602525856,
 -0.40156211504605935826451863042788609622,   0.086551159192319548209296799380484667575,   0.052343251870673974514245787611602525856,
 0.086551159192319548209296799380484667575,   -0.40156211504605935826451863042788609622,   0.052343251870673974514245787611602525856,
 -0.40156211504605935826451863042788609622,   -0.68498904414626018994477816895259857135,   0.052343251870673974514245787611602525856,
 -0.68498904414626018994477816895259857135,   -0.40156211504605935826451863042788609622,   0.052343251870673974514245787611602525856}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 18 simplex
    tempNumIPs = 67;
    nIPMap[cubDataKey(2, 18, simplex)] = tempNumIPs;
    // nIP 67
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,   0.072711470602853332335988402976151320416,
 -0.20008874386484753663926775646038722347,   -0.59982251227030492672146448707922555305,   0.066608940066780270570831020513249877437,
 -0.59982251227030492672146448707922555305,   -0.20008874386484753663926775646038722347,   0.066608940066780270570831020513249877437,
 -0.20008874386484753663926775646038722347,   -0.20008874386484753663926775646038722347,   0.066608940066780270570831020513249877437,
-0.024839396850260875373628934599793544657,   -0.95032120629947824925274213080041291069,   0.024093295267999420236853581306803869378,
 -0.95032120629947824925274213080041291069,  -0.024839396850260875373628934599793544657,   0.024093295267999420236853581306803869378,
-0.024839396850260875373628934599793544657,  -0.024839396850260875373628934599793544657,   0.024093295267999420236853581306803869378,
-0.076380987187101537262574124416325864243,   -0.84723802562579692547485175116734827151,   0.037898343013557732052184724503582760819,
 -0.84723802562579692547485175116734827151,  -0.076380987187101537262574124416325864243,   0.037898343013557732052184724503582760819,
-0.076380987187101537262574124416325864243,  -0.076380987187101537262574124416325864243,   0.037898343013557732052184724503582760819,
 -0.51547059497145607559012153340572104878,   0.030941189942912151180243066811442097565,   0.072950178817887274068289295917294690431,
 0.030941189942912151180243066811442097565,   -0.51547059497145607559012153340572104878,   0.072950178817887274068289295917294690431,
 -0.51547059497145607559012153340572104878,   -0.51547059497145607559012153340572104878,   0.072950178817887274068289295917294690431,
 -0.92233948782262881057955835940173820882,    0.84467897564525762115911671880347641764,   0.014258652039437941355528309610938991063,
  0.84467897564525762115911671880347641764,   -0.92233948782262881057955835940173820882,   0.014258652039437941355528309610938991063,
 -0.92233948782262881057955835940173820882,   -0.92233948782262881057955835940173820882,   0.014258652039437941355528309610938991063,
 -0.81610451575671360652823121099663693088,    0.63220903151342721305646242199327386177,   0.033118319904006498336116111979261819383,
  0.63220903151342721305646242199327386177,   -0.81610451575671360652823121099663693088,   0.033118319904006498336116111979261819383,
 -0.81610451575671360652823121099663693088,   -0.81610451575671360652823121099663693088,   0.033118319904006498336116111979261819383,
 -0.90839016828027843774281581498808452091,    0.54074475242935042645861991501168424548,    0.02751923246988441047583904641836977339,
  0.54074475242935042645861991501168424548,   -0.90839016828027843774281581498808452091,    0.02751923246988441047583904641836977339,
 -0.63235458414907198871580410002359972457,    0.54074475242935042645861991501168424548,    0.02751923246988441047583904641836977339,
  0.54074475242935042645861991501168424548,   -0.63235458414907198871580410002359972457,    0.02751923246988441047583904641836977339,
 -0.63235458414907198871580410002359972457,   -0.90839016828027843774281581498808452091,    0.02751923246988441047583904641836977339,
 -0.90839016828027843774281581498808452091,   -0.63235458414907198871580410002359972457,    0.02751923246988441047583904641836977339,
 -0.58730148513232413316130962964884514978,     0.3419079703884690263223113955340021876,   0.047563821800305659903183394223175954568,
   0.3419079703884690263223113955340021876,   -0.58730148513232413316130962964884514978,   0.047563821800305659903183394223175954568,
 -0.75460648525614489316100176588515703782,     0.3419079703884690263223113955340021876,   0.047563821800305659903183394223175954568,
   0.3419079703884690263223113955340021876,   -0.75460648525614489316100176588515703782,   0.047563821800305659903183394223175954568,
 -0.75460648525614489316100176588515703782,   -0.58730148513232413316130962964884514978,   0.047563821800305659903183394223175954568,
 -0.58730148513232413316130962964884514978,   -0.75460648525614489316100176588515703782,   0.047563821800305659903183394223175954568,
 -0.99220477793305323492615661731437867596,    0.20083790926851382953172593624080189463,  0.0090610690045141307984091089591295584041,
  0.20083790926851382953172593624080189463,   -0.99220477793305323492615661731437867596,  0.0090610690045141307984091089591295584041,
 -0.20863313133546059460556931892642321867,    0.20083790926851382953172593624080189463,  0.0090610690045141307984091089591295584041,
  0.20083790926851382953172593624080189463,   -0.20863313133546059460556931892642321867,  0.0090610690045141307984091089591295584041,
 -0.20863313133546059460556931892642321867,   -0.99220477793305323492615661731437867596,  0.0090610690045141307984091089591295584041,
 -0.99220477793305323492615661731437867596,   -0.20863313133546059460556931892642321867,  0.0090610690045141307984091089591295584041,
 -0.97307596651711002157127978980742349944,     0.7566843789350434355117007033556566479,   0.013680220239214362980455516076649575226,
   0.7566843789350434355117007033556566479,   -0.97307596651711002157127978980742349944,   0.013680220239214362980455516076649575226,
 -0.78360841241793341394042091354823314847,     0.7566843789350434355117007033556566479,   0.013680220239214362980455516076649575226,
   0.7566843789350434355117007033556566479,   -0.78360841241793341394042091354823314847,   0.013680220239214362980455516076649575226,
 -0.78360841241793341394042091354823314847,   -0.97307596651711002157127978980742349944,   0.013680220239214362980455516076649575226,
 -0.97307596651711002157127978980742349944,   -0.78360841241793341394042091354823314847,   0.013680220239214362980455516076649575226,
 -0.91947943306018387371123279519490851215,    0.27997618400942918982209258667637955233,   0.035494978204040810232789924461690967346,
  0.27997618400942918982209258667637955233,   -0.91947943306018387371123279519490851215,   0.035494978204040810232789924461690967346,
 -0.36049675094924531611085979148147104018,    0.27997618400942918982209258667637955233,   0.035494978204040810232789924461690967346,
  0.27997618400942918982209258667637955233,   -0.36049675094924531611085979148147104018,   0.035494978204040810232789924461690967346,
 -0.36049675094924531611085979148147104018,   -0.91947943306018387371123279519490851215,   0.035494978204040810232789924461690967346,
 -0.91947943306018387371123279519490851215,   -0.36049675094924531611085979148147104018,   0.035494978204040810232789924461690967346,
 -0.98940332962678046939842579820752941279,    0.51785895971039698809029890766872251942,   0.010021321749159443407741271393882052956,
  0.51785895971039698809029890766872251942,   -0.98940332962678046939842579820752941279,   0.010021321749159443407741271393882052956,
 -0.52845563008361651869187310946119310664,    0.51785895971039698809029890766872251942,   0.010021321749159443407741271393882052956,
  0.51785895971039698809029890766872251942,   -0.52845563008361651869187310946119310664,   0.010021321749159443407741271393882052956,
 -0.52845563008361651869187310946119310664,   -0.98940332962678046939842579820752941279,   0.010021321749159443407741271393882052956,
 -0.98940332962678046939842579820752941279,   -0.52845563008361651869187310946119310664,   0.010021321749159443407741271393882052956,
 -0.99890327991591536201516398045062601334,    0.94472145792559133356359726011634553554,  0.0024458962539221795584360481949957423674,
  0.94472145792559133356359726011634553554,   -0.99890327991591536201516398045062601334,  0.0024458962539221795584360481949957423674,
  -0.9458181780096759715484332796657195222,    0.94472145792559133356359726011634553554,  0.0024458962539221795584360481949957423674,
  0.94472145792559133356359726011634553554,    -0.9458181780096759715484332796657195222,  0.0024458962539221795584360481949957423674,
  -0.9458181780096759715484332796657195222,   -0.99890327991591536201516398045062601334,  0.0024458962539221795584360481949957423674,
 -0.99890327991591536201516398045062601334,    -0.9458181780096759715484332796657195222,  0.0024458962539221795584360481949957423674,
 -0.75882460967215071364642989357871249414,   0.091837550772389199617388042151854332785,   0.050964350623648878943912767860515151418,
 0.091837550772389199617388042151854332785,   -0.75882460967215071364642989357871249414,   0.050964350623648878943912767860515151418,
 -0.33301294110023848597095814857314183865,   0.091837550772389199617388042151854332785,   0.050964350623648878943912767860515151418,
 0.091837550772389199617388042151854332785,   -0.33301294110023848597095814857314183865,   0.050964350623648878943912767860515151418,
 -0.33301294110023848597095814857314183865,   -0.75882460967215071364642989357871249414,   0.050964350623648878943912767860515151418,
 -0.75882460967215071364642989357871249414,   -0.33301294110023848597095814857314183865,   0.050964350623648878943912767860515151418}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 18 and 19 orthotope
    tempNumIPs = 72;
    nIPMap[cubDataKey(2, 18, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 19, orthotope)] = tempNumIPs;
    // nIP 72
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0.71434425957279423719800751930789820023,                                           0,   0.097737364692828748243518350785527543045,
                                         0,    0.71434425957279423719800751930789820023,   0.097737364692828748243518350785527543045,
 -0.71434425957279423719800751930789820023,                                           0,   0.097737364692828748243518350785527543045,
                                         0,   -0.71434425957279423719800751930789820023,   0.097737364692828748243518350785527543045,
  0.26567205212096375634521276949089468244,                                           0,    0.13930761292248224570420924010933858068,
                                         0,    0.26567205212096375634521276949089468244,    0.13930761292248224570420924010933858068,
 -0.26567205212096375634521276949089468244,                                           0,    0.13930761292248224570420924010933858068,
                                         0,   -0.26567205212096375634521276949089468244,    0.13930761292248224570420924010933858068,
   0.9644342692579672514758983345396097952,                                           0,   0.034869583491887910002162414690733291135,
                                         0,     0.9644342692579672514758983345396097952,   0.034869583491887910002162414690733291135,
  -0.9644342692579672514758983345396097952,                                           0,   0.034869583491887910002162414690733291135,
                                         0,    -0.9644342692579672514758983345396097952,   0.034869583491887910002162414690733291135,
  0.80337972946768500184569945317555062404,    0.80337972946768500184569945317555062404,   0.047804547162795794795004833111469267529,
  0.80337972946768500184569945317555062404,   -0.80337972946768500184569945317555062404,   0.047804547162795794795004833111469267529,
 -0.80337972946768500184569945317555062404,    0.80337972946768500184569945317555062404,   0.047804547162795794795004833111469267529,
 -0.80337972946768500184569945317555062404,   -0.80337972946768500184569945317555062404,   0.047804547162795794795004833111469267529,
  0.99216540589713472189744527407633022081,    0.99216540589713472189744527407633022081,  0.0016177151779117608109927450603284034943,
  0.99216540589713472189744527407633022081,   -0.99216540589713472189744527407633022081,  0.0016177151779117608109927450603284034943,
 -0.99216540589713472189744527407633022081,    0.99216540589713472189744527407633022081,  0.0016177151779117608109927450603284034943,
 -0.99216540589713472189744527407633022081,   -0.99216540589713472189744527407633022081,  0.0016177151779117608109927450603284034943,
  0.92944960279960938592453972561789815525,    0.92944960279960938592453972561789815525,   0.017441048034355762464098424359798214084,
  0.92944960279960938592453972561789815525,   -0.92944960279960938592453972561789815525,   0.017441048034355762464098424359798214084,
 -0.92944960279960938592453972561789815525,    0.92944960279960938592453972561789815525,   0.017441048034355762464098424359798214084,
 -0.92944960279960938592453972561789815525,   -0.92944960279960938592453972561789815525,   0.017441048034355762464098424359798214084,
  0.51027825736934335058280520856084337461,    0.26664031459456217566653623598727229551,    0.11772583284005608002278985852480623731,
  0.26664031459456217566653623598727229551,    0.51027825736934335058280520856084337461,    0.11772583284005608002278985852480623731,
  0.51027825736934335058280520856084337461,   -0.26664031459456217566653623598727229551,    0.11772583284005608002278985852480623731,
 -0.26664031459456217566653623598727229551,    0.51027825736934335058280520856084337461,    0.11772583284005608002278985852480623731,
 -0.51027825736934335058280520856084337461,    0.26664031459456217566653623598727229551,    0.11772583284005608002278985852480623731,
  0.26664031459456217566653623598727229551,   -0.51027825736934335058280520856084337461,    0.11772583284005608002278985852480623731,
 -0.51027825736934335058280520856084337461,   -0.26664031459456217566653623598727229551,    0.11772583284005608002278985852480623731,
 -0.26664031459456217566653623598727229551,   -0.51027825736934335058280520856084337461,    0.11772583284005608002278985852480623731,
  0.39073420577524975124346417118562879486,    0.98789293314173530957079377090965768822,   0.016179577611657849789658841248953982496,
  0.98789293314173530957079377090965768822,    0.39073420577524975124346417118562879486,   0.016179577611657849789658841248953982496,
  0.39073420577524975124346417118562879486,   -0.98789293314173530957079377090965768822,   0.016179577611657849789658841248953982496,
 -0.98789293314173530957079377090965768822,    0.39073420577524975124346417118562879486,   0.016179577611657849789658841248953982496,
 -0.39073420577524975124346417118562879486,    0.98789293314173530957079377090965768822,   0.016179577611657849789658841248953982496,
  0.98789293314173530957079377090965768822,   -0.39073420577524975124346417118562879486,   0.016179577611657849789658841248953982496,
 -0.39073420577524975124346417118562879486,   -0.98789293314173530957079377090965768822,   0.016179577611657849789658841248953982496,
 -0.98789293314173530957079377090965768822,   -0.39073420577524975124346417118562879486,   0.016179577611657849789658841248953982496,
  0.71716792130974519606308556126612093997,    0.51249187721609765575479280145710018395,   0.082846618983404893940813472711113790502,
  0.51249187721609765575479280145710018395,    0.71716792130974519606308556126612093997,   0.082846618983404893940813472711113790502,
  0.71716792130974519606308556126612093997,   -0.51249187721609765575479280145710018395,   0.082846618983404893940813472711113790502,
 -0.51249187721609765575479280145710018395,    0.71716792130974519606308556126612093997,   0.082846618983404893940813472711113790502,
 -0.71716792130974519606308556126612093997,    0.51249187721609765575479280145710018395,   0.082846618983404893940813472711113790502,
  0.51249187721609765575479280145710018395,   -0.71716792130974519606308556126612093997,   0.082846618983404893940813472711113790502,
 -0.71716792130974519606308556126612093997,   -0.51249187721609765575479280145710018395,   0.082846618983404893940813472711113790502,
 -0.51249187721609765575479280145710018395,   -0.71716792130974519606308556126612093997,   0.082846618983404893940813472711113790502,
   0.2654400078112959960719970586153501497,    0.86890243415450423858119610210847838329,   0.064989081492596046079059800766815612157,
  0.86890243415450423858119610210847838329,     0.2654400078112959960719970586153501497,   0.064989081492596046079059800766815612157,
   0.2654400078112959960719970586153501497,   -0.86890243415450423858119610210847838329,   0.064989081492596046079059800766815612157,
 -0.86890243415450423858119610210847838329,     0.2654400078112959960719970586153501497,   0.064989081492596046079059800766815612157,
  -0.2654400078112959960719970586153501497,    0.86890243415450423858119610210847838329,   0.064989081492596046079059800766815612157,
  0.86890243415450423858119610210847838329,    -0.2654400078112959960719970586153501497,   0.064989081492596046079059800766815612157,
  -0.2654400078112959960719970586153501497,   -0.86890243415450423858119610210847838329,   0.064989081492596046079059800766815612157,
 -0.86890243415450423858119610210847838329,    -0.2654400078112959960719970586153501497,   0.064989081492596046079059800766815612157,
  0.62003539869325639760371263809929928256,    0.92630295580712924244748529919286127317,   0.038980552396410534776396621703868595797,
  0.92630295580712924244748529919286127317,    0.62003539869325639760371263809929928256,   0.038980552396410534776396621703868595797,
  0.62003539869325639760371263809929928256,   -0.92630295580712924244748529919286127317,   0.038980552396410534776396621703868595797,
 -0.92630295580712924244748529919286127317,    0.62003539869325639760371263809929928256,   0.038980552396410534776396621703868595797,
 -0.62003539869325639760371263809929928256,    0.92630295580712924244748529919286127317,   0.038980552396410534776396621703868595797,
  0.92630295580712924244748529919286127317,   -0.62003539869325639760371263809929928256,   0.038980552396410534776396621703868595797,
 -0.62003539869325639760371263809929928256,   -0.92630295580712924244748529919286127317,   0.038980552396410534776396621703868595797,
 -0.92630295580712924244748529919286127317,   -0.62003539869325639760371263809929928256,   0.038980552396410534776396621703868595797,
  0.80167158471859684259236358940451923912,    0.98844653068397375168746156789890756035,  0.0098894009347434843812884009858441317495,
  0.98844653068397375168746156789890756035,    0.80167158471859684259236358940451923912,  0.0098894009347434843812884009858441317495,
  0.80167158471859684259236358940451923912,   -0.98844653068397375168746156789890756035,  0.0098894009347434843812884009858441317495,
 -0.98844653068397375168746156789890756035,    0.80167158471859684259236358940451923912,  0.0098894009347434843812884009858441317495,
 -0.80167158471859684259236358940451923912,    0.98844653068397375168746156789890756035,  0.0098894009347434843812884009858441317495,
  0.98844653068397375168746156789890756035,   -0.80167158471859684259236358940451923912,  0.0098894009347434843812884009858441317495,
 -0.80167158471859684259236358940451923912,   -0.98844653068397375168746156789890756035,  0.0098894009347434843812884009858441317495,
 -0.98844653068397375168746156789890756035,   -0.80167158471859684259236358940451923912,  0.0098894009347434843812884009858441317495}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 19 simplex
    tempNumIPs = 73;
    nIPMap[cubDataKey(2, 19, simplex)] = tempNumIPs;
    // nIP 73
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,   0.068938795408024669812026416913714150021,
 -0.89495221929758206597048698795512662293,    0.78990443859516413194097397591025324585,   0.014218513195596261075438194944381222518,
  0.78990443859516413194097397591025324585,   -0.89495221929758206597048698795512662293,   0.014218513195596261075438194944381222518,
 -0.89495221929758206597048698795512662293,   -0.89495221929758206597048698795512662293,   0.014218513195596261075438194944381222518,
-0.014974649917326261006204672820750524368,   -0.97005070016534747798759065435849895126,   0.020643510285888563379375271511299695768,
 -0.97005070016534747798759065435849895126,  -0.014974649917326261006204672820750524368,   0.020643510285888563379375271511299695768,
-0.014974649917326261006204672820750524368,  -0.014974649917326261006204672820750524368,   0.020643510285888563379375271511299695768,
 -0.77710225335395724709458940859057135024,    0.55420450670791449418917881718114270048,   0.030468702186036599162431617970615364518,
  0.55420450670791449418917881718114270048,   -0.77710225335395724709458940859057135024,   0.030468702186036599162431617970615364518,
 -0.77710225335395724709458940859057135024,   -0.77710225335395724709458940859057135024,   0.030468702186036599162431617970615364518,
 -0.08161159792091268505961214209837404408,   -0.83677680415817462988077571580325191184,   0.045967180053483212460620337545594040128,
 -0.83677680415817462988077571580325191184,   -0.08161159792091268505961214209837404408,   0.045967180053483212460620337545594040128,
 -0.08161159792091268505961214209837404408,   -0.08161159792091268505961214209837404408,   0.045967180053483212460620337545594040128,
 -0.19206055489619760678392498292764093084,   -0.61587889020760478643215003414471813831,   0.063075069786309950153050389763287229858,
 -0.61587889020760478643215003414471813831,   -0.19206055489619760678392498292764093084,   0.063075069786309950153050389763287229858,
 -0.19206055489619760678392498292764093084,   -0.19206055489619760678392498292764093084,   0.063075069786309950153050389763287229858,
 -0.64365979043647137831445631799854042799,    0.28731958087294275662891263599708085598,    0.04930382969638170658102716032808733157,
  0.28731958087294275662891263599708085598,   -0.64365979043647137831445631799854042799,    0.04930382969638170658102716032808733157,
 -0.64365979043647137831445631799854042799,   -0.64365979043647137831445631799854042799,    0.04930382969638170658102716032808733157,
 -0.97672107763242110615802660011074523649,    0.95344215526484221231605320022149047297,  0.0035306455528856953194485177898948363749,
  0.95344215526484221231605320022149047297,   -0.97672107763242110615802660011074523649,  0.0035306455528856953194485177898948363749,
 -0.97672107763242110615802660011074523649,   -0.97672107763242110615802660011074523649,  0.0035306455528856953194485177898948363749,
 -0.48967673417278460498250110195791340279,  -0.020646531654430790034997796084173194424,   0.063506038732006153941897444999340725061,
-0.020646531654430790034997796084173194424,   -0.48967673417278460498250110195791340279,   0.063506038732006153941897444999340725061,
 -0.48967673417278460498250110195791340279,   -0.48967673417278460498250110195791340279,   0.063506038732006153941897444999340725061,
 -0.73860464746393519086078211530574473521,     0.6603129288005507482114383730999627094,   0.019390968973710091597859466457133882581,
   0.6603129288005507482114383730999627094,   -0.73860464746393519086078211530574473521,   0.019390968973710091597859466457133882581,
 -0.92170828133661555735065625779421797419,     0.6603129288005507482114383730999627094,   0.019390968973710091597859466457133882581,
   0.6603129288005507482114383730999627094,   -0.92170828133661555735065625779421797419,   0.019390968973710091597859466457133882581,
 -0.92170828133661555735065625779421797419,   -0.73860464746393519086078211530574473521,   0.019390968973710091597859466457133882581,
 -0.73860464746393519086078211530574473521,   -0.92170828133661555735065625779421797419,   0.019390968973710091597859466457133882581,
 -0.37736474038091749444973870042603089225,    0.11873961144060182831770759671302796842,   0.052692643954781466370362997690064596362,
  0.11873961144060182831770759671302796842,   -0.37736474038091749444973870042603089225,   0.052692643954781466370362997690064596362,
 -0.74137487105968433386796889628699707617,    0.11873961144060182831770759671302796842,   0.052692643954781466370362997690064596362,
  0.11873961144060182831770759671302796842,   -0.74137487105968433386796889628699707617,   0.052692643954781466370362997690064596362,
 -0.74137487105968433386796889628699707617,   -0.37736474038091749444973870042603089225,   0.052692643954781466370362997690064596362,
 -0.37736474038091749444973870042603089225,   -0.74137487105968433386796889628699707617,   0.052692643954781466370362997690064596362,
 -0.99586214820679038477978859628118904207,    0.26662658625756826609373531798200113825,   0.006564153103671638511890683462332859992,
  0.26662658625756826609373531798200113825,   -0.99586214820679038477978859628118904207,   0.006564153103671638511890683462332859992,
 -0.27076443805077788131394672170081209618,    0.26662658625756826609373531798200113825,   0.006564153103671638511890683462332859992,
  0.26662658625756826609373531798200113825,   -0.27076443805077788131394672170081209618,   0.006564153103671638511890683462332859992,
 -0.27076443805077788131394672170081209618,   -0.99586214820679038477978859628118904207,   0.006564153103671638511890683462332859992,
 -0.99586214820679038477978859628118904207,   -0.27076443805077788131394672170081209618,   0.006564153103671638511890683462332859992,
 -0.85087941079674664292586974917469643186,    0.40800963993208421976175838156331845374,   0.036215889862424900129305127397133399259,
  0.40800963993208421976175838156331845374,   -0.85087941079674664292586974917469643186,   0.036215889862424900129305127397133399259,
 -0.55713022913533757683588863238862202188,    0.40800963993208421976175838156331845374,   0.036215889862424900129305127397133399259,
  0.40800963993208421976175838156331845374,   -0.55713022913533757683588863238862202188,   0.036215889862424900129305127397133399259,
 -0.55713022913533757683588863238862202188,   -0.85087941079674664292586974917469643186,   0.036215889862424900129305127397133399259,
 -0.85087941079674664292586974917469643186,   -0.55713022913533757683588863238862202188,   0.036215889862424900129305127397133399259,
 -0.98998542348529101773301781123483422136,    0.70513390875377839935650663034800986154,  0.0058526302069404008122936576981139743993,
  0.70513390875377839935650663034800986154,   -0.98998542348529101773301781123483422136,  0.0058526302069404008122936576981139743993,
 -0.71514848526848738162348881911317564017,    0.70513390875377839935650663034800986154,  0.0058526302069404008122936576981139743993,
  0.70513390875377839935650663034800986154,   -0.71514848526848738162348881911317564017,  0.0058526302069404008122936576981139743993,
 -0.71514848526848738162348881911317564017,   -0.98998542348529101773301781123483422136,  0.0058526302069404008122936576981139743993,
 -0.98998542348529101773301781123483422136,   -0.71514848526848738162348881911317564017,  0.0058526302069404008122936576981139743993,
 -0.91822397760796624772174420894241076741,    0.21016795813741592822569883837786491478,   0.032204325528048216419984026965287582219,
  0.21016795813741592822569883837786491478,   -0.91822397760796624772174420894241076741,   0.032204325528048216419984026965287582219,
 -0.29194398052944968050395462943545414738,    0.21016795813741592822569883837786491478,   0.032204325528048216419984026965287582219,
  0.21016795813741592822569883837786491478,   -0.29194398052944968050395462943545414738,   0.032204325528048216419984026965287582219,
 -0.29194398052944968050395462943545414738,   -0.91822397760796624772174420894241076741,   0.032204325528048216419984026965287582219,
 -0.91822397760796624772174420894241076741,   -0.29194398052944968050395462943545414738,   0.032204325528048216419984026965287582219,
 -0.51621084207884086284394200561371723534,    0.48636273791487271949314123346299601262,   0.016911774999072993635042140907181887112,
  0.48636273791487271949314123346299601262,   -0.51621084207884086284394200561371723534,   0.016911774999072993635042140907181887112,
 -0.97015189583603185664919922784927877729,    0.48636273791487271949314123346299601262,   0.016911774999072993635042140907181887112,
  0.48636273791487271949314123346299601262,   -0.97015189583603185664919922784927877729,   0.016911774999072993635042140907181887112,
 -0.97015189583603185664919922784927877729,   -0.51621084207884086284394200561371723534,   0.016911774999072993635042140907181887112,
 -0.51621084207884086284394200561371723534,   -0.97015189583603185664919922784927877729,   0.016911774999072993635042140907181887112,
 -0.87982744935538660401520750808428413826,    0.86027539775361029331257605015603810029,  0.0066544027257187765179463625108825701742,
  0.86027539775361029331257605015603810029,   -0.87982744935538660401520750808428413826,  0.0066544027257187765179463625108825701742,
 -0.98044794839822368929736854207175396203,    0.86027539775361029331257605015603810029,  0.0066544027257187765179463625108825701742,
  0.86027539775361029331257605015603810029,   -0.98044794839822368929736854207175396203,  0.0066544027257187765179463625108825701742,
 -0.98044794839822368929736854207175396203,   -0.87982744935538660401520750808428413826,  0.0066544027257187765179463625108825701742,
 -0.87982744935538660401520750808428413826,   -0.98044794839822368929736854207175396203,  0.0066544027257187765179463625108825701742}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 20 simplex
    tempNumIPs = 79;
    nIPMap[cubDataKey(2, 20, simplex)] = tempNumIPs;
    // nIP 79
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,   0.055640442805812463631566788253898051532,
 -0.49084146465332177523029470374523072315,  -0.018317070693356449539410592509538553699,    0.05633280523008099012864178412473069936,
-0.018317070693356449539410592509538553699,   -0.49084146465332177523029470374523072315,    0.05633280523008099012864178412473069936,
 -0.49084146465332177523029470374523072315,   -0.49084146465332177523029470374523072315,    0.05633280523008099012864178412473069936,
 -0.97804771794320447192070397628224037507,    0.95609543588640894384140795256448075015,  0.0031953631642664794478945847330657770722,
  0.95609543588640894384140795256448075015,   -0.97804771794320447192070397628224037507,  0.0031953631642664794478945847330657770722,
 -0.97804771794320447192070397628224037507,   -0.97804771794320447192070397628224037507,  0.0031953631642664794478945847330657770722,
 -0.78123280657657080596625505020033275589,    0.56246561315314161193251010040066551178,   0.031320923104298133684407003648894203447,
  0.56246561315314161193251010040066551178,   -0.78123280657657080596625505020033275589,   0.031320923104298133684407003648894203447,
 -0.78123280657657080596625505020033275589,   -0.78123280657657080596625505020033275589,   0.031320923104298133684407003648894203447,
 -0.62741000451091811445653498388475129105,     0.2548200090218362289130699677695025821,   0.036693851897011657919252854118138043537,
   0.2548200090218362289130699677695025821,   -0.62741000451091811445653498388475129105,   0.036693851897011657919252854118138043537,
 -0.62741000451091811445653498388475129105,   -0.62741000451091811445653498388475129105,   0.036693851897011657919252854118138043537,
 -0.10889788608815036962555544313157118467,   -0.78220422782369926074888911373685763067,   0.037809599732929790747314162640799596721,
 -0.78220422782369926074888911373685763067,   -0.10889788608815036962555544313157118467,   0.037809599732929790747314162640799596721,
 -0.10889788608815036962555544313157118467,   -0.10889788608815036962555544313157118467,   0.037809599732929790747314162640799596721,
 -0.92537823880223061201299256742556755612,    0.85075647760446122402598513485113511224,   0.008645101642662310104692722157806774508,
  0.85075647760446122402598513485113511224,   -0.92537823880223061201299256742556755612,   0.008645101642662310104692722157806774508,
 -0.92537823880223061201299256742556755612,   -0.92537823880223061201299256742556755612,   0.008645101642662310104692722157806774508,
 -0.21314930436580028264207251299344481704,   -0.57370139126839943471585497401311036593,   0.055152202516281836053550890863277393079,
 -0.57370139126839943471585497401311036593,   -0.21314930436580028264207251299344481704,   0.055152202516281836053550890863277393079,
 -0.21314930436580028264207251299344481704,   -0.21314930436580028264207251299344481704,   0.055152202516281836053550890863277393079,
-0.047508776919001973579296255534818088061,   -0.90498244616199605284140748893036382388,   0.028407301213633761967680465924479290382,
 -0.90498244616199605284140748893036382388,  -0.047508776919001973579296255534818088061,   0.028407301213633761967680465924479290382,
-0.047508776919001973579296255534818088061,  -0.047508776919001973579296255534818088061,   0.028407301213633761967680465924479290382,
 -0.98485843899060694347353308294287690802,    0.66659102367647249434430238008617319241,  0.0088115896742339902561312708157853846482,
  0.66659102367647249434430238008617319241,   -0.98485843899060694347353308294287690802,  0.0088115896742339902561312708157853846482,
 -0.68173258468586555087076929714329628438,    0.66659102367647249434430238008617319241,  0.0088115896742339902561312708157853846482,
  0.66659102367647249434430238008617319241,   -0.68173258468586555087076929714329628438,  0.0088115896742339902561312708157853846482,
 -0.68173258468586555087076929714329628438,   -0.98485843899060694347353308294287690802,  0.0088115896742339902561312708157853846482,
 -0.98485843899060694347353308294287690802,   -0.68173258468586555087076929714329628438,  0.0088115896742339902561312708157853846482,
 -0.90687927018467136697246089078946388707,     0.5098430057270950102124174847737658903,   0.023945594315818760080585256144999165493,
   0.5098430057270950102124174847737658903,   -0.90687927018467136697246089078946388707,   0.023945594315818760080585256144999165493,
 -0.60296373554242364323995659398430200323,     0.5098430057270950102124174847737658903,   0.023945594315818760080585256144999165493,
   0.5098430057270950102124174847737658903,   -0.60296373554242364323995659398430200323,   0.023945594315818760080585256144999165493,
 -0.60296373554242364323995659398430200323,   -0.90687927018467136697246089078946388707,   0.023945594315818760080585256144999165493,
 -0.90687927018467136697246089078946388707,   -0.60296373554242364323995659398430200323,   0.023945594315818760080585256144999165493,
 -0.87181882878313187990026866476433962457,    0.86210895356788437332353836160651207109,  0.0045194784085034621179235150583929560485,
  0.86210895356788437332353836160651207109,   -0.87181882878313187990026866476433962457,  0.0045194784085034621179235150583929560485,
 -0.99029012478475249342326969684217244652,    0.86210895356788437332353836160651207109,  0.0045194784085034621179235150583929560485,
  0.86210895356788437332353836160651207109,   -0.99029012478475249342326969684217244652,  0.0045194784085034621179235150583929560485,
 -0.99029012478475249342326969684217244652,   -0.87181882878313187990026866476433962457,  0.0045194784085034621179235150583929560485,
 -0.87181882878313187990026866476433962457,   -0.99029012478475249342326969684217244652,  0.0045194784085034621179235150583929560485,
 -0.89002504171402637873154550296837045792,    0.22375540709485139987629530764328543032,   0.034668902268877332248730626058035734221,
  0.22375540709485139987629530764328543032,   -0.89002504171402637873154550296837045792,   0.034668902268877332248730626058035734221,
  -0.3337303653808250211447498046749149724,    0.22375540709485139987629530764328543032,   0.034668902268877332248730626058035734221,
  0.22375540709485139987629530764328543032,    -0.3337303653808250211447498046749149724,   0.034668902268877332248730626058035734221,
  -0.3337303653808250211447498046749149724,   -0.89002504171402637873154550296837045792,   0.034668902268877332248730626058035734221,
 -0.89002504171402637873154550296837045792,    -0.3337303653808250211447498046749149724,   0.034668902268877332248730626058035734221,
  -0.8000954074237226892341843533057676038,    0.72336803787297349848896164841146460112,   0.016582846110455430851049450635630476598,
  0.72336803787297349848896164841146460112,    -0.8000954074237226892341843533057676038,   0.016582846110455430851049450635630476598,
 -0.92327263044925080925477729510569699731,    0.72336803787297349848896164841146460112,   0.016582846110455430851049450635630476598,
  0.72336803787297349848896164841146460112,   -0.92327263044925080925477729510569699731,   0.016582846110455430851049450635630476598,
 -0.92327263044925080925477729510569699731,    -0.8000954074237226892341843533057676038,   0.016582846110455430851049450635630476598,
  -0.8000954074237226892341843533057676038,   -0.92327263044925080925477729510569699731,   0.016582846110455430851049450635630476598,
 -0.78754559055945991527489139589221510388,    0.35633147577927110921599545999123193124,   0.030890431288396919377481650650853417183,
  0.35633147577927110921599545999123193124,   -0.78754559055945991527489139589221510388,   0.030890431288396919377481650650853417183,
 -0.56878588521981119394110406409901682736,    0.35633147577927110921599545999123193124,   0.030890431288396919377481650650853417183,
  0.35633147577927110921599545999123193124,   -0.56878588521981119394110406409901682736,   0.030890431288396919377481650650853417183,
 -0.56878588521981119394110406409901682736,   -0.78754559055945991527489139589221510388,   0.030890431288396919377481650650853417183,
 -0.78754559055945991527489139589221510388,   -0.56878588521981119394110406409901682736,   0.030890431288396919377481650650853417183,
  -0.1599524823675518407071141661054621128,     0.1402893857819467195500916946480886286,    0.01478272600102119191287687445277038169,
   0.1402893857819467195500916946480886286,    -0.1599524823675518407071141661054621128,    0.01478272600102119191287687445277038169,
  -0.9803369034143948788429775285426265158,     0.1402893857819467195500916946480886286,    0.01478272600102119191287687445277038169,
   0.1402893857819467195500916946480886286,    -0.9803369034143948788429775285426265158,    0.01478272600102119191287687445277038169,
  -0.9803369034143948788429775285426265158,    -0.1599524823675518407071141661054621128,    0.01478272600102119191287687445277038169,
  -0.1599524823675518407071141661054621128,    -0.9803369034143948788429775285426265158,    0.01478272600102119191287687445277038169,
 -0.36427975232845596804806465673948070356,   0.084663608344856170859761688802646336641,   0.046766982927310947730184310937781849814,
 0.084663608344856170859761688802646336641,   -0.36427975232845596804806465673948070356,   0.046766982927310947730184310937781849814,
 -0.72038385601640020281169703206316563308,   0.084663608344856170859761688802646336641,   0.046766982927310947730184310937781849814,
 0.084663608344856170859761688802646336641,   -0.72038385601640020281169703206316563308,   0.046766982927310947730184310937781849814,
 -0.72038385601640020281169703206316563308,   -0.36427975232845596804806465673948070356,   0.046766982927310947730184310937781849814,
 -0.36427975232845596804806465673948070356,   -0.72038385601640020281169703206316563308,   0.046766982927310947730184310937781849814,
 -0.97852557428797782533895087798420844581,    0.41736275144064735566397781689608203229,   0.014312800953830741459725346431171736663,
  0.41736275144064735566397781689608203229,   -0.97852557428797782533895087798420844581,   0.014312800953830741459725346431171736663,
 -0.43883717715266953032502693891187358648,    0.41736275144064735566397781689608203229,   0.014312800953830741459725346431171736663,
  0.41736275144064735566397781689608203229,   -0.43883717715266953032502693891187358648,   0.014312800953830741459725346431171736663,
 -0.43883717715266953032502693891187358648,   -0.97852557428797782533895087798420844581,   0.014312800953830741459725346431171736663,
 -0.97852557428797782533895087798420844581,   -0.43883717715266953032502693891187358648,   0.014312800953830741459725346431171736663}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 20 and 21 orthotope
    tempNumIPs = 85;
    nIPMap[cubDataKey(2, 20, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(2, 21, orthotope)] = tempNumIPs;
    // nIP 85
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0,                                           0,    0.13519284035614282527922840385434746149,
  0.47335107435822415202690901838844666219,                                           0,     0.1067941587859403825037542148829566842,
                                         0,    0.47335107435822415202690901838844666219,     0.1067941587859403825037542148829566842,
 -0.47335107435822415202690901838844666219,                                           0,     0.1067941587859403825037542148829566842,
                                         0,   -0.47335107435822415202690901838844666219,     0.1067941587859403825037542148829566842,
  0.83527842975586830160520893009495292054,                                           0,   0.066257208004944385159031059020122440339,
                                         0,    0.83527842975586830160520893009495292054,   0.066257208004944385159031059020122440339,
 -0.83527842975586830160520893009495292054,                                           0,   0.066257208004944385159031059020122440339,
                                         0,   -0.83527842975586830160520893009495292054,   0.066257208004944385159031059020122440339,
   0.2573719006072289819053755213317636919,     0.2573719006072289819053755213317636919,    0.11988448364634046970472217835662908028,
   0.2573719006072289819053755213317636919,    -0.2573719006072289819053755213317636919,    0.11988448364634046970472217835662908028,
  -0.2573719006072289819053755213317636919,     0.2573719006072289819053755213317636919,    0.11988448364634046970472217835662908028,
  -0.2573719006072289819053755213317636919,    -0.2573719006072289819053755213317636919,    0.11988448364634046970472217835662908028,
  0.96333783211562335381484352624945730322,    0.96333783211562335381484352624945730322,  0.0082294612892200088282217117254895783997,
  0.96333783211562335381484352624945730322,   -0.96333783211562335381484352624945730322,  0.0082294612892200088282217117254895783997,
 -0.96333783211562335381484352624945730322,    0.96333783211562335381484352624945730322,  0.0082294612892200088282217117254895783997,
 -0.96333783211562335381484352624945730322,   -0.96333783211562335381484352624945730322,  0.0082294612892200088282217117254895783997,
  0.86245192537965155140267607631968215083,    0.86245192537965155140267607631968215083,   0.030603640928775651038141155704768806683,
  0.86245192537965155140267607631968215083,   -0.86245192537965155140267607631968215083,   0.030603640928775651038141155704768806683,
 -0.86245192537965155140267607631968215083,    0.86245192537965155140267607631968215083,   0.030603640928775651038141155704768806683,
 -0.86245192537965155140267607631968215083,   -0.86245192537965155140267607631968215083,   0.030603640928775651038141155704768806683,
  0.49689796251934570364793617005129236026,    0.49689796251934570364793617005129236026,   0.096791791893595210842497711127833756611,
  0.49689796251934570364793617005129236026,   -0.49689796251934570364793617005129236026,   0.096791791893595210842497711127833756611,
 -0.49689796251934570364793617005129236026,    0.49689796251934570364793617005129236026,   0.096791791893595210842497711127833756611,
 -0.49689796251934570364793617005129236026,   -0.49689796251934570364793617005129236026,   0.096791791893595210842497711127833756611,
  0.70433217519540056955932591780080720749,    0.70433217519540056955932591780080720749,   0.061516341481316561201749944009961586508,
  0.70433217519540056955932591780080720749,   -0.70433217519540056955932591780080720749,   0.061516341481316561201749944009961586508,
 -0.70433217519540056955932591780080720749,    0.70433217519540056955932591780080720749,   0.061516341481316561201749944009961586508,
 -0.70433217519540056955932591780080720749,   -0.70433217519540056955932591780080720749,   0.061516341481316561201749944009961586508,
  0.24187818547670197670997054949192246819,    0.67414621995531773950245618191533516236,   0.086728499513208224903936696999552926752,
  0.67414621995531773950245618191533516236,    0.24187818547670197670997054949192246819,   0.086728499513208224903936696999552926752,
  0.24187818547670197670997054949192246819,   -0.67414621995531773950245618191533516236,   0.086728499513208224903936696999552926752,
 -0.67414621995531773950245618191533516236,    0.24187818547670197670997054949192246819,   0.086728499513208224903936696999552926752,
 -0.24187818547670197670997054949192246819,    0.67414621995531773950245618191533516236,   0.086728499513208224903936696999552926752,
  0.67414621995531773950245618191533516236,   -0.24187818547670197670997054949192246819,   0.086728499513208224903936696999552926752,
 -0.24187818547670197670997054949192246819,   -0.67414621995531773950245618191533516236,   0.086728499513208224903936696999552926752,
 -0.67414621995531773950245618191533516236,   -0.24187818547670197670997054949192246819,   0.086728499513208224903936696999552926752,
  0.48015696631279513715726389107913072158,    0.82464737527092073827881221191700729562,   0.055879117404127348063047858204117827292,
  0.82464737527092073827881221191700729562,    0.48015696631279513715726389107913072158,   0.055879117404127348063047858204117827292,
  0.48015696631279513715726389107913072158,   -0.82464737527092073827881221191700729562,   0.055879117404127348063047858204117827292,
 -0.82464737527092073827881221191700729562,    0.48015696631279513715726389107913072158,   0.055879117404127348063047858204117827292,
 -0.48015696631279513715726389107913072158,    0.82464737527092073827881221191700729562,   0.055879117404127348063047858204117827292,
  0.82464737527092073827881221191700729562,   -0.48015696631279513715726389107913072158,   0.055879117404127348063047858204117827292,
 -0.48015696631279513715726389107913072158,   -0.82464737527092073827881221191700729562,   0.055879117404127348063047858204117827292,
 -0.82464737527092073827881221191700729562,   -0.48015696631279513715726389107913072158,   0.055879117404127348063047858204117827292,
  0.93224193592175396207500211941993080129,    0.27069918410166491089101469370373344757,   0.037128827440913817653330718954876960924,
  0.27069918410166491089101469370373344757,    0.93224193592175396207500211941993080129,   0.037128827440913817653330718954876960924,
  0.93224193592175396207500211941993080129,   -0.27069918410166491089101469370373344757,   0.037128827440913817653330718954876960924,
 -0.27069918410166491089101469370373344757,    0.93224193592175396207500211941993080129,   0.037128827440913817653330718954876960924,
 -0.93224193592175396207500211941993080129,    0.27069918410166491089101469370373344757,   0.037128827440913817653330718954876960924,
  0.27069918410166491089101469370373344757,   -0.93224193592175396207500211941993080129,   0.037128827440913817653330718954876960924,
 -0.93224193592175396207500211941993080129,   -0.27069918410166491089101469370373344757,   0.037128827440913817653330718954876960924,
 -0.27069918410166491089101469370373344757,   -0.93224193592175396207500211941993080129,   0.037128827440913817653330718954876960924,
  0.93490232582401054613354438740771124115,    0.67503956743707526937208896355984125883,   0.028775600851020531767579367626940765984,
  0.67503956743707526937208896355984125883,    0.93490232582401054613354438740771124115,   0.028775600851020531767579367626940765984,
  0.93490232582401054613354438740771124115,   -0.67503956743707526937208896355984125883,   0.028775600851020531767579367626940765984,
 -0.67503956743707526937208896355984125883,    0.93490232582401054613354438740771124115,   0.028775600851020531767579367626940765984,
 -0.93490232582401054613354438740771124115,    0.67503956743707526937208896355984125883,   0.028775600851020531767579367626940765984,
  0.67503956743707526937208896355984125883,   -0.93490232582401054613354438740771124115,   0.028775600851020531767579367626940765984,
 -0.93490232582401054613354438740771124115,   -0.67503956743707526937208896355984125883,   0.028775600851020531767579367626940765984,
 -0.67503956743707526937208896355984125883,   -0.93490232582401054613354438740771124115,   0.028775600851020531767579367626940765984,
  0.99045546752952415272197458895882111757,    0.48891625117716687937811887553934444444,   0.011509520442493167578261138308188460848,
  0.48891625117716687937811887553934444444,    0.99045546752952415272197458895882111757,   0.011509520442493167578261138308188460848,
  0.99045546752952415272197458895882111757,   -0.48891625117716687937811887553934444444,   0.011509520442493167578261138308188460848,
 -0.48891625117716687937811887553934444444,    0.99045546752952415272197458895882111757,   0.011509520442493167578261138308188460848,
 -0.99045546752952415272197458895882111757,    0.48891625117716687937811887553934444444,   0.011509520442493167578261138308188460848,
  0.48891625117716687937811887553934444444,   -0.99045546752952415272197458895882111757,   0.011509520442493167578261138308188460848,
 -0.99045546752952415272197458895882111757,   -0.48891625117716687937811887553934444444,   0.011509520442493167578261138308188460848,
 -0.48891625117716687937811887553934444444,   -0.99045546752952415272197458895882111757,   0.011509520442493167578261138308188460848,
  0.83113524597590138585974389212898407353,    0.98896061138657239166243578378631145048,  0.0075284821304016457091249399363317448912,
  0.98896061138657239166243578378631145048,    0.83113524597590138585974389212898407353,  0.0075284821304016457091249399363317448912,
  0.83113524597590138585974389212898407353,   -0.98896061138657239166243578378631145048,  0.0075284821304016457091249399363317448912,
 -0.98896061138657239166243578378631145048,    0.83113524597590138585974389212898407353,  0.0075284821304016457091249399363317448912,
 -0.83113524597590138585974389212898407353,    0.98896061138657239166243578378631145048,  0.0075284821304016457091249399363317448912,
  0.98896061138657239166243578378631145048,   -0.83113524597590138585974389212898407353,  0.0075284821304016457091249399363317448912,
 -0.83113524597590138585974389212898407353,   -0.98896061138657239166243578378631145048,  0.0075284821304016457091249399363317448912,
 -0.98896061138657239166243578378631145048,   -0.83113524597590138585974389212898407353,  0.0075284821304016457091249399363317448912,
 0.091469600922291296668443828573471880948,    0.98401714848959892866013755872564589893,   0.010512304158251076525756742074316914114,
  0.98401714848959892866013755872564589893,   0.091469600922291296668443828573471880948,   0.010512304158251076525756742074316914114,
 0.091469600922291296668443828573471880948,   -0.98401714848959892866013755872564589893,   0.010512304158251076525756742074316914114,
 -0.98401714848959892866013755872564589893,   0.091469600922291296668443828573471880948,   0.010512304158251076525756742074316914114,
-0.091469600922291296668443828573471880948,    0.98401714848959892866013755872564589893,   0.010512304158251076525756742074316914114,
  0.98401714848959892866013755872564589893,  -0.091469600922291296668443828573471880948,   0.010512304158251076525756742074316914114,
-0.091469600922291296668443828573471880948,   -0.98401714848959892866013755872564589893,   0.010512304158251076525756742074316914114,
 -0.98401714848959892866013755872564589893,  -0.091469600922291296668443828573471880948,   0.010512304158251076525756742074316914114}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //dimension 3
    tempDim = 3;
    //order 1 orthotope
    tempNumIPs = 1;
    nIPMap[cubDataKey(tempDim, 1, orthotope)] = tempNumIPs;
    // nIP 1
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {0,  0,  0,  8}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 1 simplex
    tempNumIPs = 1;
    nIPMap[cubDataKey(tempDim, 1, simplex)] = tempNumIPs;
    // nIP 1
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.5,                                     -0.5,                                     -0.5,  1.3333333333333333333333333333333333333}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 2 simplex
    tempNumIPs = 4;
    nIPMap[cubDataKey(tempDim, 2, simplex)] = tempNumIPs;
    // nIP 4
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.72360679774997896964091736687312762354,  -0.72360679774997896964091736687312762354,   0.17082039324993690892275210061938287063,   0.33333333333333333333333333333333333333,
-0.72360679774997896964091736687312762354,   0.17082039324993690892275210061938287063,  -0.72360679774997896964091736687312762354,   0.33333333333333333333333333333333333333,
 0.17082039324993690892275210061938287063,  -0.72360679774997896964091736687312762354,  -0.72360679774997896964091736687312762354,   0.33333333333333333333333333333333333333,
-0.72360679774997896964091736687312762354,  -0.72360679774997896964091736687312762354,  -0.72360679774997896964091736687312762354,   0.33333333333333333333333333333333333333}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 2 and 3 orthotope
    tempNumIPs = 6;
    nIPMap[cubDataKey(tempDim, 2, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(tempDim, 3, orthotope)] = tempNumIPs;
    // nIP 6
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-1,                                        0,                                        0,  1.3333333333333333333333333333333333333,
                                      0,                                        0,                                        1,  1.3333333333333333333333333333333333333,
                                      0,                                        1,                                        0,  1.3333333333333333333333333333333333333,
                                      0,                                        0,                                       -1,  1.3333333333333333333333333333333333333,
                                      1,                                        0,                                        0,  1.3333333333333333333333333333333333333,
                                      0,                                       -1,                                        0,  1.3333333333333333333333333333333333333}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 3 simplex
    tempNumIPs = 8;
    nIPMap[cubDataKey(tempDim, 3, simplex)] = tempNumIPs;
    // nIP 8
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,   -0.9689798150982901207378151874892027072,   0.18162379004944980942342872025562069427,
-0.34367339496723662642072827083693243093,   -0.9689798150982901207378151874892027072,  -0.34367339496723662642072827083693243093,   0.18162379004944980942342872025562069427,
 -0.9689798150982901207378151874892027072,  -0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,   0.18162379004944980942342872025562069427,
-0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,   0.18162379004944980942342872025562069427,
-0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,   0.35171650060942837529461966476512015033,   0.15170954328388352390990461307771263906,
-0.78390550020314279176487322158837338344,   0.35171650060942837529461966476512015033,  -0.78390550020314279176487322158837338344,   0.15170954328388352390990461307771263906,
 0.35171650060942837529461966476512015033,  -0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,   0.15170954328388352390990461307771263906,
-0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,   0.15170954328388352390990461307771263906}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 4 and 5 orthotope
    tempNumIPs = 14;
    nIPMap[cubDataKey(tempDim, 4, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(tempDim, 5, orthotope)] = tempNumIPs;
    // nIP 14
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.79582242575422146326454882047613584616,                                          0,                                          0,   0.88642659279778393351800554016620498615,
                                        0,                                          0,   0.79582242575422146326454882047613584616,   0.88642659279778393351800554016620498615,
                                        0,   0.79582242575422146326454882047613584616,                                          0,   0.88642659279778393351800554016620498615,
                                        0,                                          0,  -0.79582242575422146326454882047613584616,   0.88642659279778393351800554016620498615,
 0.79582242575422146326454882047613584616,                                          0,                                          0,   0.88642659279778393351800554016620498615,
                                        0,  -0.79582242575422146326454882047613584616,                                          0,   0.88642659279778393351800554016620498615,
 0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039,
-0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039,
-0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039,
-0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039,
-0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039,
 0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039,
 0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039,
 0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 4 and 5 simplex
    tempNumIPs = 14;
    nIPMap[cubDataKey(tempDim, 4, simplex)] = tempNumIPs;
    nIPMap[cubDataKey(tempDim, 5, simplex)] = tempNumIPs;
    // nIP 14
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.37822816147339878040530853247308433401,   -0.37822816147339878040530853247308433401,   -0.86531551557980365878407440258074699796,    0.15025056762402113439891420311104844508,
 -0.37822816147339878040530853247308433401,   -0.86531551557980365878407440258074699796,   -0.37822816147339878040530853247308433401,    0.15025056762402113439891420311104844508,
 -0.86531551557980365878407440258074699796,   -0.37822816147339878040530853247308433401,   -0.37822816147339878040530853247308433401,    0.15025056762402113439891420311104844508,
 -0.37822816147339878040530853247308433401,   -0.37822816147339878040530853247308433401,   -0.37822816147339878040530853247308433401,    0.15025056762402113439891420311104844508,
 -0.81452949937821754719535217252593878951,   -0.81452949937821754719535217252593878951,    0.44358849813465264158605651757781636853,   0.097990724155149266058280273981770004697,
 -0.81452949937821754719535217252593878951,    0.44358849813465264158605651757781636853,   -0.81452949937821754719535217252593878951,   0.097990724155149266058280273981770004697,
  0.44358849813465264158605651757781636853,   -0.81452949937821754719535217252593878951,   -0.81452949937821754719535217252593878951,   0.097990724155149266058280273981770004697,
 -0.81452949937821754719535217252593878951,   -0.81452949937821754719535217252593878951,   -0.81452949937821754719535217252593878951,   0.097990724155149266058280273981770004697,
 -0.90899259174870070101623894744132112187,  -0.091007408251299298983761052558678878131,  -0.091007408251299298983761052558678878131,    0.05672802770277528858409257082700992237,
-0.091007408251299298983761052558678878131,   -0.90899259174870070101623894744132112187,  -0.091007408251299298983761052558678878131,    0.05672802770277528858409257082700992237,
 -0.90899259174870070101623894744132112187,   -0.90899259174870070101623894744132112187,  -0.091007408251299298983761052558678878131,    0.05672802770277528858409257082700992237,
 -0.90899259174870070101623894744132112187,  -0.091007408251299298983761052558678878131,   -0.90899259174870070101623894744132112187,    0.05672802770277528858409257082700992237,
-0.091007408251299298983761052558678878131,   -0.90899259174870070101623894744132112187,   -0.90899259174870070101623894744132112187,    0.05672802770277528858409257082700992237,
-0.091007408251299298983761052558678878131,  -0.091007408251299298983761052558678878131,   -0.90899259174870070101623894744132112187,    0.05672802770277528858409257082700992237}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 6 simplex
    tempNumIPs = 24;
    nIPMap[cubDataKey(tempDim, 6, simplex)] = tempNumIPs;
    // nIP 24
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.91865208293077729376884110208717988144,  -0.91865208293077729376884110208717988144,   0.75595624879233188130652330626153964433,  0.013436281407094190597350983261249151023,
-0.91865208293077729376884110208717988144,   0.75595624879233188130652330626153964433,  -0.91865208293077729376884110208717988144,  0.013436281407094190597350983261249151023,
 0.75595624879233188130652330626153964433,  -0.91865208293077729376884110208717988144,  -0.91865208293077729376884110208717988144,  0.013436281407094190597350983261249151023,
-0.91865208293077729376884110208717988144,  -0.91865208293077729376884110208717988144,  -0.91865208293077729376884110208717988144,  0.013436281407094190597350983261249151023,
-0.35532421971544897931201105847501574951,  -0.35532421971544897931201105847501574951,  -0.93402734085365306206396682457495275146,  0.073809575391539629460204370471634688549,
-0.35532421971544897931201105847501574951,  -0.93402734085365306206396682457495275146,  -0.35532421971544897931201105847501574951,  0.073809575391539629460204370471634688549,
-0.93402734085365306206396682457495275146,  -0.35532421971544897931201105847501574951,  -0.35532421971544897931201105847501574951,  0.073809575391539629460204370471634688549,
-0.35532421971544897931201105847501574951,  -0.35532421971544897931201105847501574951,  -0.35532421971544897931201105847501574951,  0.073809575391539629460204370471634688549,
 -0.5707942574816959414223215612274300172,   -0.5707942574816959414223215612274300172,  -0.28761722755491217573303531631770994839,  0.053230333677556656132920836743306636618,
 -0.5707942574816959414223215612274300172,  -0.28761722755491217573303531631770994839,   -0.5707942574816959414223215612274300172,  0.053230333677556656132920836743306636618,
-0.28761722755491217573303531631770994839,   -0.5707942574816959414223215612274300172,   -0.5707942574816959414223215612274300172,  0.053230333677556656132920836743306636618,
 -0.5707942574816959414223215612274300172,   -0.5707942574816959414223215612274300172,   -0.5707942574816959414223215612274300172,  0.053230333677556656132920836743306636618,
 0.20601132958329828273486227812187937257,  -0.87267799624996494940152894478854603924,  -0.46065533708336838393180438854478729409,  0.064285714285714285714285714285714285714,
 0.20601132958329828273486227812187937257,  -0.87267799624996494940152894478854603924,  -0.87267799624996494940152894478854603924,  0.064285714285714285714285714285714285714,
-0.87267799624996494940152894478854603924,  -0.87267799624996494940152894478854603924,   0.20601132958329828273486227812187937257,  0.064285714285714285714285714285714285714,
-0.46065533708336838393180438854478729409,   0.20601132958329828273486227812187937257,  -0.87267799624996494940152894478854603924,  0.064285714285714285714285714285714285714,
-0.87267799624996494940152894478854603924,  -0.46065533708336838393180438854478729409,   0.20601132958329828273486227812187937257,  0.064285714285714285714285714285714285714,
-0.87267799624996494940152894478854603924,   0.20601132958329828273486227812187937257,  -0.87267799624996494940152894478854603924,  0.064285714285714285714285714285714285714,
-0.46065533708336838393180438854478729409,  -0.87267799624996494940152894478854603924,   0.20601132958329828273486227812187937257,  0.064285714285714285714285714285714285714,
-0.87267799624996494940152894478854603924,  -0.46065533708336838393180438854478729409,  -0.87267799624996494940152894478854603924,  0.064285714285714285714285714285714285714,
-0.87267799624996494940152894478854603924,  -0.87267799624996494940152894478854603924,  -0.46065533708336838393180438854478729409,  0.064285714285714285714285714285714285714,
-0.87267799624996494940152894478854603924,   0.20601132958329828273486227812187937257,  -0.46065533708336838393180438854478729409,  0.064285714285714285714285714285714285714,
-0.46065533708336838393180438854478729409,  -0.87267799624996494940152894478854603924,  -0.87267799624996494940152894478854603924,  0.064285714285714285714285714285714285714,
 0.20601132958329828273486227812187937257,  -0.46065533708336838393180438854478729409,  -0.87267799624996494940152894478854603924,  0.064285714285714285714285714285714285714}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    // order 6 and 7 orthotope
    tempNumIPs = 34;
    nIPMap[cubDataKey(tempDim, 6, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(tempDim, 7, orthotope)] = tempNumIPs;
    // nIP 34
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.98808616118867654860459729594655654482,                                          0,                                          0,   0.20012989098356700695575850638756794678,
                                        0,                                          0,   0.98808616118867654860459729594655654482,   0.20012989098356700695575850638756794678,
                                        0,   0.98808616118867654860459729594655654482,                                          0,   0.20012989098356700695575850638756794678,
                                        0,                                          0,  -0.98808616118867654860459729594655654482,   0.20012989098356700695575850638756794678,
 0.98808616118867654860459729594655654482,                                          0,                                          0,   0.20012989098356700695575850638756794678,
                                        0,  -0.98808616118867654860459729594655654482,                                          0,   0.20012989098356700695575850638756794678,
 0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
-0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
-0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
-0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
-0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
 0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
 0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
 0.40795516735831936633080916251520586856,  -0.40795516735831936633080916251520586856,   0.40795516735831936633080916251520586856,   0.45715385606985274296912190785204854361,
 0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
-0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
-0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
-0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
-0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
 0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
 0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
 0.78110282100411850665136814715928740203,  -0.78110282100411850665136814715928740203,   0.78110282100411850665136814715928740203,   0.15379614006595869408326783499527494637,
-0.84805227568403872779622085759407839383,  -0.84805227568403872779622085759407839383,                                          0,   0.15930172375100887182052758490800036662,
 0.84805227568403872779622085759407839383,                                          0,  -0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662,
                                        0,   0.84805227568403872779622085759407839383,  -0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662,
 0.84805227568403872779622085759407839383,   0.84805227568403872779622085759407839383,                                          0,   0.15930172375100887182052758490800036662,
 0.84805227568403872779622085759407839383,                                          0,   0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662,
                                        0,  -0.84805227568403872779622085759407839383,   0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662,
                                        0,  -0.84805227568403872779622085759407839383,  -0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662,
-0.84805227568403872779622085759407839383,                                          0,   0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662,
-0.84805227568403872779622085759407839383,   0.84805227568403872779622085759407839383,                                          0,   0.15930172375100887182052758490800036662,
 0.84805227568403872779622085759407839383,  -0.84805227568403872779622085759407839383,                                          0,   0.15930172375100887182052758490800036662,
                                        0,   0.84805227568403872779622085759407839383,   0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662,
-0.84805227568403872779622085759407839383,                                          0,  -0.84805227568403872779622085759407839383,   0.15930172375100887182052758490800036662}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 7 simplex
    tempNumIPs = 35;
    nIPMap[cubDataKey(tempDim, 7, simplex)] = tempNumIPs;
    // nIP 35
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.5,                                       -0.5,                                       -0.5,   0.12731371928550779848077124815630184245,
-0.36859770044359440115314000081337702315,  -0.36859770044359440115314000081337702315,  -0.89420689866921679654057999755986893056,   0.05643944161328937210171489439806231665,
-0.36859770044359440115314000081337702315,  -0.89420689866921679654057999755986893056,  -0.36859770044359440115314000081337702315,   0.05643944161328937210171489439806231665,
-0.89420689866921679654057999755986893056,  -0.36859770044359440115314000081337702315,  -0.36859770044359440115314000081337702315,   0.05643944161328937210171489439806231665,
-0.36859770044359440115314000081337702315,  -0.36859770044359440115314000081337702315,  -0.36859770044359440115314000081337702315,   0.05643944161328937210171489439806231665,
-0.89902035480320726247389235402687506805,  -0.10097964519679273752610764597312493195,  -0.10097964519679273752610764597312493195,   0.04252923711047677324569976544392327613,
-0.10097964519679273752610764597312493195,  -0.89902035480320726247389235402687506805,  -0.10097964519679273752610764597312493195,   0.04252923711047677324569976544392327613,
-0.89902035480320726247389235402687506805,  -0.89902035480320726247389235402687506805,  -0.10097964519679273752610764597312493195,   0.04252923711047677324569976544392327613,
-0.89902035480320726247389235402687506805,  -0.10097964519679273752610764597312493195,  -0.89902035480320726247389235402687506805,   0.04252923711047677324569976544392327613,
-0.10097964519679273752610764597312493195,  -0.89902035480320726247389235402687506805,  -0.89902035480320726247389235402687506805,   0.04252923711047677324569976544392327613,
-0.10097964519679273752610764597312493195,  -0.10097964519679273752610764597312493195,  -0.89902035480320726247389235402687506805,   0.04252923711047677324569976544392327613,
 0.15034327517400004696648315404461503933,  -0.62233233794799790452713779229082848999,  -0.90567859927800423791220756946295805935,  0.049609507637779495159487414921974827082,
 0.15034327517400004696648315404461503933,  -0.62233233794799790452713779229082848999,  -0.62233233794799790452713779229082848999,  0.049609507637779495159487414921974827082,
-0.62233233794799790452713779229082848999,  -0.62233233794799790452713779229082848999,   0.15034327517400004696648315404461503933,  0.049609507637779495159487414921974827082,
-0.90567859927800423791220756946295805935,   0.15034327517400004696648315404461503933,  -0.62233233794799790452713779229082848999,  0.049609507637779495159487414921974827082,
-0.62233233794799790452713779229082848999,  -0.90567859927800423791220756946295805935,   0.15034327517400004696648315404461503933,  0.049609507637779495159487414921974827082,
-0.62233233794799790452713779229082848999,   0.15034327517400004696648315404461503933,  -0.62233233794799790452713779229082848999,  0.049609507637779495159487414921974827082,
-0.90567859927800423791220756946295805935,  -0.62233233794799790452713779229082848999,   0.15034327517400004696648315404461503933,  0.049609507637779495159487414921974827082,
-0.62233233794799790452713779229082848999,  -0.90567859927800423791220756946295805935,  -0.62233233794799790452713779229082848999,  0.049609507637779495159487414921974827082,
-0.62233233794799790452713779229082848999,  -0.62233233794799790452713779229082848999,  -0.90567859927800423791220756946295805935,  0.049609507637779495159487414921974827082,
-0.62233233794799790452713779229082848999,   0.15034327517400004696648315404461503933,  -0.90567859927800423791220756946295805935,  0.049609507637779495159487414921974827082,
-0.90567859927800423791220756946295805935,  -0.62233233794799790452713779229082848999,  -0.62233233794799790452713779229082848999,  0.049609507637779495159487414921974827082,
 0.15034327517400004696648315404461503933,  -0.90567859927800423791220756946295805935,  -0.62233233794799790452713779229082848999,  0.049609507637779495159487414921974827082,
 0.62166048219709712223621075969646478759,  -0.95746905491703350802232779700036011785,  -0.70672237236303010619155516569574455189,  0.010814361106537788754804577988128720209,
 0.62166048219709712223621075969646478759,  -0.95746905491703350802232779700036011785,  -0.95746905491703350802232779700036011785,  0.010814361106537788754804577988128720209,
-0.95746905491703350802232779700036011785,  -0.95746905491703350802232779700036011785,   0.62166048219709712223621075969646478759,  0.010814361106537788754804577988128720209,
-0.70672237236303010619155516569574455189,   0.62166048219709712223621075969646478759,  -0.95746905491703350802232779700036011785,  0.010814361106537788754804577988128720209,
-0.95746905491703350802232779700036011785,  -0.70672237236303010619155516569574455189,   0.62166048219709712223621075969646478759,  0.010814361106537788754804577988128720209,
-0.95746905491703350802232779700036011785,   0.62166048219709712223621075969646478759,  -0.95746905491703350802232779700036011785,  0.010814361106537788754804577988128720209,
-0.70672237236303010619155516569574455189,  -0.95746905491703350802232779700036011785,   0.62166048219709712223621075969646478759,  0.010814361106537788754804577988128720209,
-0.95746905491703350802232779700036011785,  -0.70672237236303010619155516569574455189,  -0.95746905491703350802232779700036011785,  0.010814361106537788754804577988128720209,
-0.95746905491703350802232779700036011785,  -0.95746905491703350802232779700036011785,  -0.70672237236303010619155516569574455189,  0.010814361106537788754804577988128720209,
-0.95746905491703350802232779700036011785,   0.62166048219709712223621075969646478759,  -0.70672237236303010619155516569574455189,  0.010814361106537788754804577988128720209,
-0.70672237236303010619155516569574455189,  -0.95746905491703350802232779700036011785,  -0.95746905491703350802232779700036011785,  0.010814361106537788754804577988128720209,
 0.62166048219709712223621075969646478759,  -0.70672237236303010619155516569574455189,  -0.95746905491703350802232779700036011785,  0.010814361106537788754804577988128720209}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 8 simplex
    tempNumIPs = 46;
    nIPMap[cubDataKey(tempDim, 8, simplex)] = tempNumIPs;
    // nIP 46
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.78409455007557830303561343854585276851,   -0.78409455007557830303561343854585276851,    0.35228365022673490910684031563755830553,   0.035235534544545108539097313565645827902,
 -0.78409455007557830303561343854585276851,    0.35228365022673490910684031563755830553,   -0.78409455007557830303561343854585276851,   0.035235534544545108539097313565645827902,
  0.35228365022673490910684031563755830553,   -0.78409455007557830303561343854585276851,   -0.78409455007557830303561343854585276851,   0.035235534544545108539097313565645827902,
 -0.78409455007557830303561343854585276851,   -0.78409455007557830303561343854585276851,   -0.78409455007557830303561343854585276851,   0.035235534544545108539097313565645827902,
   -0.629781024434826859446183999496503797,     -0.629781024434826859446183999496503797,     -0.110656926695519421661448001510488609,   0.069375663418318043619861410028472677606,
   -0.629781024434826859446183999496503797,     -0.110656926695519421661448001510488609,     -0.629781024434826859446183999496503797,   0.069375663418318043619861410028472677606,
   -0.110656926695519421661448001510488609,     -0.629781024434826859446183999496503797,     -0.629781024434826859446183999496503797,   0.069375663418318043619861410028472677606,
   -0.629781024434826859446183999496503797,     -0.629781024434826859446183999496503797,     -0.629781024434826859446183999496503797,   0.069375663418318043619861410028472677606,
 -0.91536691263046543675183591431855057189,   -0.91536691263046543675183591431855057189,    0.74610073789139631025550774295565171568,   0.010033674871386933040967405110292766211,
 -0.91536691263046543675183591431855057189,    0.74610073789139631025550774295565171568,   -0.91536691263046543675183591431855057189,   0.010033674871386933040967405110292766211,
  0.74610073789139631025550774295565171568,   -0.91536691263046543675183591431855057189,   -0.91536691263046543675183591431855057189,   0.010033674871386933040967405110292766211,
 -0.91536691263046543675183591431855057189,   -0.91536691263046543675183591431855057189,   -0.91536691263046543675183591431855057189,   0.010033674871386933040967405110292766211,
 -0.37163658175192200927334245312373662702,   -0.37163658175192200927334245312373662702,   -0.88509025474423397217997264062879011895,   0.055685043809246530275971813171748834223,
 -0.37163658175192200927334245312373662702,   -0.88509025474423397217997264062879011895,   -0.37163658175192200927334245312373662702,   0.055685043809246530275971813171748834223,
 -0.88509025474423397217997264062879011895,   -0.37163658175192200927334245312373662702,   -0.37163658175192200927334245312373662702,   0.055685043809246530275971813171748834223,
 -0.37163658175192200927334245312373662702,   -0.37163658175192200927334245312373662702,   -0.37163658175192200927334245312373662702,   0.055685043809246530275971813171748834223,
 -0.12881734283233958954373367640680593526,   -0.87118265716766041045626632359319406474,   -0.87118265716766041045626632359319406474,   0.048374573681745097334959362864457783412,
 -0.87118265716766041045626632359319406474,   -0.12881734283233958954373367640680593526,   -0.87118265716766041045626632359319406474,   0.048374573681745097334959362864457783412,
 -0.12881734283233958954373367640680593526,   -0.12881734283233958954373367640680593526,   -0.87118265716766041045626632359319406474,   0.048374573681745097334959362864457783412,
 -0.12881734283233958954373367640680593526,   -0.87118265716766041045626632359319406474,   -0.12881734283233958954373367640680593526,   0.048374573681745097334959362864457783412,
 -0.87118265716766041045626632359319406474,   -0.12881734283233958954373367640680593526,   -0.12881734283233958954373367640680593526,   0.048374573681745097334959362864457783412,
 -0.87118265716766041045626632359319406474,   -0.87118265716766041045626632359319406474,   -0.12881734283233958954373367640680593526,   0.048374573681745097334959362864457783412,
  0.43492812685261664658005711425116290088,   -0.95713213974573885031100182115353118963,   -0.52066384736113894595805347194410052161,  0.0095425371877925775336250150191510735196,
  0.43492812685261664658005711425116290088,   -0.95713213974573885031100182115353118963,   -0.95713213974573885031100182115353118963,  0.0095425371877925775336250150191510735196,
 -0.95713213974573885031100182115353118963,   -0.95713213974573885031100182115353118963,    0.43492812685261664658005711425116290088,  0.0095425371877925775336250150191510735196,
 -0.52066384736113894595805347194410052161,    0.43492812685261664658005711425116290088,   -0.95713213974573885031100182115353118963,  0.0095425371877925775336250150191510735196,
 -0.95713213974573885031100182115353118963,   -0.52066384736113894595805347194410052161,    0.43492812685261664658005711425116290088,  0.0095425371877925775336250150191510735196,
 -0.95713213974573885031100182115353118963,    0.43492812685261664658005711425116290088,   -0.95713213974573885031100182115353118963,  0.0095425371877925775336250150191510735196,
 -0.52066384736113894595805347194410052161,   -0.95713213974573885031100182115353118963,    0.43492812685261664658005711425116290088,  0.0095425371877925775336250150191510735196,
 -0.95713213974573885031100182115353118963,   -0.52066384736113894595805347194410052161,   -0.95713213974573885031100182115353118963,  0.0095425371877925775336250150191510735196,
 -0.95713213974573885031100182115353118963,   -0.95713213974573885031100182115353118963,   -0.52066384736113894595805347194410052161,  0.0095425371877925775336250150191510735196,
 -0.95713213974573885031100182115353118963,    0.43492812685261664658005711425116290088,   -0.52066384736113894595805347194410052161,  0.0095425371877925775336250150191510735196,
 -0.52066384736113894595805347194410052161,   -0.95713213974573885031100182115353118963,   -0.95713213974573885031100182115353118963,  0.0095425371877925775336250150191510735196,
  0.43492812685261664658005711425116290088,   -0.52066384736113894595805347194410052161,   -0.95713213974573885031100182115353118963,  0.0095425371877925775336250150191510735196,
  0.16759475660428881185333950841466005283,   -0.59172133224794175913941940243963088242,   -0.98415209210840529357450070353539828799,   0.020604648201280446418040434034344443905,
  0.16759475660428881185333950841466005283,   -0.59172133224794175913941940243963088242,   -0.59172133224794175913941940243963088242,   0.020604648201280446418040434034344443905,
 -0.59172133224794175913941940243963088242,   -0.59172133224794175913941940243963088242,    0.16759475660428881185333950841466005283,   0.020604648201280446418040434034344443905,
 -0.98415209210840529357450070353539828799,    0.16759475660428881185333950841466005283,   -0.59172133224794175913941940243963088242,   0.020604648201280446418040434034344443905,
 -0.59172133224794175913941940243963088242,   -0.98415209210840529357450070353539828799,    0.16759475660428881185333950841466005283,   0.020604648201280446418040434034344443905,
 -0.59172133224794175913941940243963088242,    0.16759475660428881185333950841466005283,   -0.59172133224794175913941940243963088242,   0.020604648201280446418040434034344443905,
 -0.98415209210840529357450070353539828799,   -0.59172133224794175913941940243963088242,    0.16759475660428881185333950841466005283,   0.020604648201280446418040434034344443905,
 -0.59172133224794175913941940243963088242,   -0.98415209210840529357450070353539828799,   -0.59172133224794175913941940243963088242,   0.020604648201280446418040434034344443905,
 -0.59172133224794175913941940243963088242,   -0.59172133224794175913941940243963088242,   -0.98415209210840529357450070353539828799,   0.020604648201280446418040434034344443905,
 -0.59172133224794175913941940243963088242,    0.16759475660428881185333950841466005283,   -0.98415209210840529357450070353539828799,   0.020604648201280446418040434034344443905,
 -0.98415209210840529357450070353539828799,   -0.59172133224794175913941940243963088242,   -0.59172133224794175913941940243963088242,   0.020604648201280446418040434034344443905,
  0.16759475660428881185333950841466005283,   -0.98415209210840529357450070353539828799,   -0.59172133224794175913941940243963088242,   0.020604648201280446418040434034344443905}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    // order 8 and 9 orthotope
    tempNumIPs = 58;
    nIPMap[cubDataKey(tempDim, 8, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(tempDim, 9, orthotope)] = tempNumIPs;
    // nIP 58
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.61368146959170899383488488974055688452,                                          0,                                          0,   0.43327499574965454300983079319432200813,
                                        0,                                          0,   0.61368146959170899383488488974055688452,   0.43327499574965454300983079319432200813,
                                        0,   0.61368146959170899383488488974055688452,                                          0,   0.43327499574965454300983079319432200813,
                                        0,                                          0,  -0.61368146959170899383488488974055688452,   0.43327499574965454300983079319432200813,
 0.61368146959170899383488488974055688452,                                          0,                                          0,   0.43327499574965454300983079319432200813,
                                        0,  -0.61368146959170899383488488974055688452,                                          0,   0.43327499574965454300983079319432200813,
 0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
-0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
-0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
-0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
-0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
 0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
 0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
 0.87009978466197591761506380886392483306,  -0.87009978466197591761506380886392483306,   0.87009978466197591761506380886392483306,  0.050148795299349029867451487724661687604,
 0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
-0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
-0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
-0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
-0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
 0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
 0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
 0.56411080702003005426661899866307283066,  -0.56411080702003005426661899866307283066,   0.56411080702003005426661899866307283066,   0.19885983814402350032086871858560920797,
-0.87768712325767828648677575899433236264,  -0.87768712325767828648677575899433236264,                                          0,  0.091789806136177642171244588919646331302,
 0.87768712325767828648677575899433236264,                                          0,  -0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
                                        0,   0.87768712325767828648677575899433236264,  -0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
 0.87768712325767828648677575899433236264,   0.87768712325767828648677575899433236264,                                          0,  0.091789806136177642171244588919646331302,
 0.87768712325767828648677575899433236264,                                          0,   0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
                                        0,  -0.87768712325767828648677575899433236264,   0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
                                        0,  -0.87768712325767828648677575899433236264,  -0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
-0.87768712325767828648677575899433236264,                                          0,   0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
-0.87768712325767828648677575899433236264,   0.87768712325767828648677575899433236264,                                          0,  0.091789806136177642171244588919646331302,
 0.87768712325767828648677575899433236264,  -0.87768712325767828648677575899433236264,                                          0,  0.091789806136177642171244588919646331302,
                                        0,   0.87768712325767828648677575899433236264,   0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
-0.87768712325767828648677575899433236264,                                          0,  -0.87768712325767828648677575899433236264,  0.091789806136177642171244588919646331302,
 0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
-0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
 0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
-0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
 0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
 0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
-0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.43226790263086216441602486151694383909,   0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
-0.43226790263086216441602486151694383909,   0.93853042186467174532897686960307879508,   0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379,
 0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,  -0.93853042186467174532897686960307879508,   0.09611680351337336643247993847150603379,
-0.93853042186467174532897686960307879508,  -0.43226790263086216441602486151694383909,  -0.43226790263086216441602486151694383909,   0.09611680351337336643247993847150603379}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
    //order 9 simplex
    tempNumIPs = 59;
    nIPMap[cubDataKey(tempDim, 9, simplex)] = tempNumIPs;
    // nIP 59
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.5,                                         -0.5,                                         -0.5,     0.07734739854997367644741906257553729918,
  -0.99999999876036601109069831529161985423,    -0.99999999876036601109069831529161985423,      0.9999999962810980332720949458748595627,  8.5759042345675186085903030258742437834e-05,
  -0.99999999876036601109069831529161985423,      0.9999999962810980332720949458748595627,    -0.99999999876036601109069831529161985423,  8.5759042345675186085903030258742437834e-05,
    0.9999999962810980332720949458748595627,    -0.99999999876036601109069831529161985423,    -0.99999999876036601109069831529161985423,  8.5759042345675186085903030258742437834e-05,
  -0.99999999876036601109069831529161985423,    -0.99999999876036601109069831529161985423,    -0.99999999876036601109069831529161985423,  8.5759042345675186085903030258742437834e-05,
  -0.67845092920947681169601712222921603853,    -0.67845092920947681169601712222921603853,    0.035352787628430435088051366687648115583,    0.030897784616567277310532826869162311937,
  -0.67845092920947681169601712222921603853,    0.035352787628430435088051366687648115583,    -0.67845092920947681169601712222921603853,    0.030897784616567277310532826869162311937,
  0.035352787628430435088051366687648115583,    -0.67845092920947681169601712222921603853,    -0.67845092920947681169601712222921603853,    0.030897784616567277310532826869162311937,
  -0.67845092920947681169601712222921603853,    -0.67845092920947681169601712222921603853,    -0.67845092920947681169601712222921603853,    0.030897784616567277310532826869162311937,
  -0.35544695635715805159186380319222134027,    -0.35544695635715805159186380319222134027,    -0.93365913092852584522440859042333597919,    0.039417216447239047523644877985518239107,
  -0.35544695635715805159186380319222134027,    -0.93365913092852584522440859042333597919,    -0.35544695635715805159186380319222134027,    0.039417216447239047523644877985518239107,
  -0.93365913092852584522440859042333597919,    -0.35544695635715805159186380319222134027,    -0.35544695635715805159186380319222134027,    0.039417216447239047523644877985518239107,
  -0.35544695635715805159186380319222134027,    -0.35544695635715805159186380319222134027,    -0.35544695635715805159186380319222134027,    0.039417216447239047523644877985518239107,
  -0.90978216330917283407807385501200614331,    -0.90978216330917283407807385501200614331,     0.72934648992751850223422156503601842994,    0.010751973306154910354337519688619045729,
  -0.90978216330917283407807385501200614331,     0.72934648992751850223422156503601842994,    -0.90978216330917283407807385501200614331,    0.010751973306154910354337519688619045729,
   0.72934648992751850223422156503601842994,    -0.90978216330917283407807385501200614331,    -0.90978216330917283407807385501200614331,    0.010751973306154910354337519688619045729,
  -0.90978216330917283407807385501200614331,    -0.90978216330917283407807385501200614331,    -0.90978216330917283407807385501200614331,    0.010751973306154910354337519688619045729,
  -0.77540690799124790934402070093968875729,    -0.22459309200875209065597929906031124271,    -0.22459309200875209065597929906031124271,     0.05084544013826995406889748131136091475,
  -0.22459309200875209065597929906031124271,    -0.77540690799124790934402070093968875729,    -0.22459309200875209065597929906031124271,     0.05084544013826995406889748131136091475,
  -0.77540690799124790934402070093968875729,    -0.77540690799124790934402070093968875729,    -0.22459309200875209065597929906031124271,     0.05084544013826995406889748131136091475,
  -0.77540690799124790934402070093968875729,    -0.22459309200875209065597929906031124271,    -0.77540690799124790934402070093968875729,     0.05084544013826995406889748131136091475,
  -0.22459309200875209065597929906031124271,    -0.77540690799124790934402070093968875729,    -0.77540690799124790934402070093968875729,     0.05084544013826995406889748131136091475,
  -0.22459309200875209065597929906031124271,    -0.22459309200875209065597929906031124271,    -0.77540690799124790934402070093968875729,     0.05084544013826995406889748131136091475,
    -0.994890841533917338064711949103255315,   -0.082257102495081453456435781598997046246,    -0.84059495347591975502241648769875059251,    0.011179229597731402927583520512290878612,
    -0.994890841533917338064711949103255315,   -0.082257102495081453456435781598997046246,   -0.082257102495081453456435781598997046246,    0.011179229597731402927583520512290878612,
 -0.082257102495081453456435781598997046246,   -0.082257102495081453456435781598997046246,      -0.994890841533917338064711949103255315,    0.011179229597731402927583520512290878612,
  -0.84059495347591975502241648769875059251,      -0.994890841533917338064711949103255315,   -0.082257102495081453456435781598997046246,    0.011179229597731402927583520512290878612,
 -0.082257102495081453456435781598997046246,    -0.84059495347591975502241648769875059251,      -0.994890841533917338064711949103255315,    0.011179229597731402927583520512290878612,
 -0.082257102495081453456435781598997046246,      -0.994890841533917338064711949103255315,   -0.082257102495081453456435781598997046246,    0.011179229597731402927583520512290878612,
  -0.84059495347591975502241648769875059251,   -0.082257102495081453456435781598997046246,      -0.994890841533917338064711949103255315,    0.011179229597731402927583520512290878612,
 -0.082257102495081453456435781598997046246,    -0.84059495347591975502241648769875059251,   -0.082257102495081453456435781598997046246,    0.011179229597731402927583520512290878612,
 -0.082257102495081453456435781598997046246,   -0.082257102495081453456435781598997046246,    -0.84059495347591975502241648769875059251,    0.011179229597731402927583520512290878612,
 -0.082257102495081453456435781598997046246,      -0.994890841533917338064711949103255315,    -0.84059495347591975502241648769875059251,    0.011179229597731402927583520512290878612,
  -0.84059495347591975502241648769875059251,   -0.082257102495081453456435781598997046246,   -0.082257102495081453456435781598997046246,    0.011179229597731402927583520512290878612,
    -0.994890841533917338064711949103255315,    -0.84059495347591975502241648769875059251,   -0.082257102495081453456435781598997046246,    0.011179229597731402927583520512290878612,
   0.43670065288414901811354448912928750012,    -0.93244825862932284418902737202159413719,    -0.57180413562550332973548974508609922574,    0.013646079136993770600501763121325612648,
   0.43670065288414901811354448912928750012,    -0.93244825862932284418902737202159413719,    -0.93244825862932284418902737202159413719,    0.013646079136993770600501763121325612648,
  -0.93244825862932284418902737202159413719,    -0.93244825862932284418902737202159413719,     0.43670065288414901811354448912928750012,    0.013646079136993770600501763121325612648,
  -0.57180413562550332973548974508609922574,     0.43670065288414901811354448912928750012,    -0.93244825862932284418902737202159413719,    0.013646079136993770600501763121325612648,
  -0.93244825862932284418902737202159413719,    -0.57180413562550332973548974508609922574,     0.43670065288414901811354448912928750012,    0.013646079136993770600501763121325612648,
  -0.93244825862932284418902737202159413719,     0.43670065288414901811354448912928750012,    -0.93244825862932284418902737202159413719,    0.013646079136993770600501763121325612648,
  -0.57180413562550332973548974508609922574,    -0.93244825862932284418902737202159413719,     0.43670065288414901811354448912928750012,    0.013646079136993770600501763121325612648,
  -0.93244825862932284418902737202159413719,    -0.57180413562550332973548974508609922574,    -0.93244825862932284418902737202159413719,    0.013646079136993770600501763121325612648,
  -0.93244825862932284418902737202159413719,    -0.93244825862932284418902737202159413719,    -0.57180413562550332973548974508609922574,    0.013646079136993770600501763121325612648,
  -0.93244825862932284418902737202159413719,     0.43670065288414901811354448912928750012,    -0.57180413562550332973548974508609922574,    0.013646079136993770600501763121325612648,
  -0.57180413562550332973548974508609922574,    -0.93244825862932284418902737202159413719,    -0.93244825862932284418902737202159413719,    0.013646079136993770600501763121325612648,
   0.43670065288414901811354448912928750012,    -0.57180413562550332973548974508609922574,    -0.93244825862932284418902737202159413719,    0.013646079136993770600501763121325612648,
  -0.93116817884364945982158957577713628669,    -0.63271726038014422042206258028726154683,     0.19660269960393790066571473635165938035,    0.027366554623984184053091789082666607808,
  -0.93116817884364945982158957577713628669,    -0.63271726038014422042206258028726154683,    -0.63271726038014422042206258028726154683,    0.027366554623984184053091789082666607808,
  -0.63271726038014422042206258028726154683,    -0.63271726038014422042206258028726154683,    -0.93116817884364945982158957577713628669,    0.027366554623984184053091789082666607808,
   0.19660269960393790066571473635165938035,    -0.93116817884364945982158957577713628669,    -0.63271726038014422042206258028726154683,    0.027366554623984184053091789082666607808,
  -0.63271726038014422042206258028726154683,     0.19660269960393790066571473635165938035,    -0.93116817884364945982158957577713628669,    0.027366554623984184053091789082666607808,
  -0.63271726038014422042206258028726154683,    -0.93116817884364945982158957577713628669,    -0.63271726038014422042206258028726154683,    0.027366554623984184053091789082666607808,
   0.19660269960393790066571473635165938035,    -0.63271726038014422042206258028726154683,    -0.93116817884364945982158957577713628669,    0.027366554623984184053091789082666607808,
  -0.63271726038014422042206258028726154683,     0.19660269960393790066571473635165938035,    -0.63271726038014422042206258028726154683,    0.027366554623984184053091789082666607808,
  -0.63271726038014422042206258028726154683,    -0.63271726038014422042206258028726154683,     0.19660269960393790066571473635165938035,    0.027366554623984184053091789082666607808,
  -0.63271726038014422042206258028726154683,    -0.93116817884364945982158957577713628669,     0.19660269960393790066571473635165938035,    0.027366554623984184053091789082666607808,
   0.19660269960393790066571473635165938035,    -0.63271726038014422042206258028726154683,    -0.63271726038014422042206258028726154683,    0.027366554623984184053091789082666607808,
  -0.93116817884364945982158957577713628669,     0.19660269960393790066571473635165938035,    -0.63271726038014422042206258028726154683,    0.027366554623984184053091789082666607808}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    //order 10 simplex
    tempNumIPs = 81;
    nIPMap[cubDataKey(tempDim, 10, simplex)] = tempNumIPs;
    // nIP 81
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.5,                                         -0.5,                                         -0.5,    0.063199698074694317846318428237401466134,
   -0.3754998626096227045403833626263450896,     -0.3754998626096227045403833626263450896,    -0.87350041217113188637884991212096473119,    0.035916079989691599737018881339842778615,
   -0.3754998626096227045403833626263450896,    -0.87350041217113188637884991212096473119,     -0.3754998626096227045403833626263450896,    0.035916079989691599737018881339842778615,
  -0.87350041217113188637884991212096473119,     -0.3754998626096227045403833626263450896,     -0.3754998626096227045403833626263450896,    0.035916079989691599737018881339842778615,
   -0.3754998626096227045403833626263450896,     -0.3754998626096227045403833626263450896,     -0.3754998626096227045403833626263450896,    0.035916079989691599737018881339842778615,
  -0.77138069228530769882525760469269910942,    -0.77138069228530769882525760469269910942,     0.31414207685592309647577281407809732827,    0.013158879622391177646076980573564101693,
  -0.77138069228530769882525760469269910942,     0.31414207685592309647577281407809732827,    -0.77138069228530769882525760469269910942,    0.013158879622391177646076980573564101693,
   0.31414207685592309647577281407809732827,    -0.77138069228530769882525760469269910942,    -0.77138069228530769882525760469269910942,    0.013158879622391177646076980573564101693,
  -0.77138069228530769882525760469269910942,    -0.77138069228530769882525760469269910942,    -0.77138069228530769882525760469269910942,    0.013158879622391177646076980573564101693,
  -0.66902794876077789679101975111094717095,    -0.17913852156206901142420431149697662502,    -0.97269500811508408036057162589509957901,    0.015191841626926975498161246507619114195,
  -0.66902794876077789679101975111094717095,    -0.17913852156206901142420431149697662502,    -0.17913852156206901142420431149697662502,    0.015191841626926975498161246507619114195,
  -0.17913852156206901142420431149697662502,    -0.17913852156206901142420431149697662502,    -0.66902794876077789679101975111094717095,    0.015191841626926975498161246507619114195,
  -0.97269500811508408036057162589509957901,    -0.66902794876077789679101975111094717095,    -0.17913852156206901142420431149697662502,    0.015191841626926975498161246507619114195,
  -0.17913852156206901142420431149697662502,    -0.97269500811508408036057162589509957901,    -0.66902794876077789679101975111094717095,    0.015191841626926975498161246507619114195,
  -0.17913852156206901142420431149697662502,    -0.66902794876077789679101975111094717095,    -0.17913852156206901142420431149697662502,    0.015191841626926975498161246507619114195,
  -0.97269500811508408036057162589509957901,    -0.17913852156206901142420431149697662502,    -0.66902794876077789679101975111094717095,    0.015191841626926975498161246507619114195,
  -0.17913852156206901142420431149697662502,    -0.97269500811508408036057162589509957901,    -0.17913852156206901142420431149697662502,    0.015191841626926975498161246507619114195,
  -0.17913852156206901142420431149697662502,    -0.17913852156206901142420431149697662502,    -0.97269500811508408036057162589509957901,    0.015191841626926975498161246507619114195,
  -0.17913852156206901142420431149697662502,    -0.66902794876077789679101975111094717095,    -0.97269500811508408036057162589509957901,    0.015191841626926975498161246507619114195,
  -0.97269500811508408036057162589509957901,    -0.17913852156206901142420431149697662502,    -0.17913852156206901142420431149697662502,    0.015191841626926975498161246507619114195,
  -0.66902794876077789679101975111094717095,    -0.97269500811508408036057162589509957901,    -0.17913852156206901142420431149697662502,    0.015191841626926975498161246507619114195,
    0.8859775346904097323952611738365015259,    -0.98772398235041850430481257350316929785,     -0.9105295699895727237856360268301629302,  0.00048259245911900483231983784641134908616,
    0.8859775346904097323952611738365015259,    -0.98772398235041850430481257350316929785,    -0.98772398235041850430481257350316929785,  0.00048259245911900483231983784641134908616,
  -0.98772398235041850430481257350316929785,    -0.98772398235041850430481257350316929785,      0.8859775346904097323952611738365015259,  0.00048259245911900483231983784641134908616,
   -0.9105295699895727237856360268301629302,      0.8859775346904097323952611738365015259,    -0.98772398235041850430481257350316929785,  0.00048259245911900483231983784641134908616,
  -0.98772398235041850430481257350316929785,     -0.9105295699895727237856360268301629302,      0.8859775346904097323952611738365015259,  0.00048259245911900483231983784641134908616,
  -0.98772398235041850430481257350316929785,      0.8859775346904097323952611738365015259,    -0.98772398235041850430481257350316929785,  0.00048259245911900483231983784641134908616,
   -0.9105295699895727237856360268301629302,    -0.98772398235041850430481257350316929785,      0.8859775346904097323952611738365015259,  0.00048259245911900483231983784641134908616,
  -0.98772398235041850430481257350316929785,     -0.9105295699895727237856360268301629302,    -0.98772398235041850430481257350316929785,  0.00048259245911900483231983784641134908616,
  -0.98772398235041850430481257350316929785,    -0.98772398235041850430481257350316929785,     -0.9105295699895727237856360268301629302,  0.00048259245911900483231983784641134908616,
  -0.98772398235041850430481257350316929785,      0.8859775346904097323952611738365015259,     -0.9105295699895727237856360268301629302,  0.00048259245911900483231983784641134908616,
   -0.9105295699895727237856360268301629302,    -0.98772398235041850430481257350316929785,    -0.98772398235041850430481257350316929785,  0.00048259245911900483231983784641134908616,
    0.8859775346904097323952611738365015259,     -0.9105295699895727237856360268301629302,    -0.98772398235041850430481257350316929785,  0.00048259245911900483231983784641134908616,
 -0.045619240191439298911787183406185558751,    -0.75789963770882114801220999680989894474,    -0.43858148439091840506379282297401655176,    0.034319642640608095038714683012872940214,
 -0.045619240191439298911787183406185558751,    -0.75789963770882114801220999680989894474,    -0.75789963770882114801220999680989894474,    0.034319642640608095038714683012872940214,
  -0.75789963770882114801220999680989894474,    -0.75789963770882114801220999680989894474,   -0.045619240191439298911787183406185558751,    0.034319642640608095038714683012872940214,
  -0.43858148439091840506379282297401655176,   -0.045619240191439298911787183406185558751,    -0.75789963770882114801220999680989894474,    0.034319642640608095038714683012872940214,
  -0.75789963770882114801220999680989894474,    -0.43858148439091840506379282297401655176,   -0.045619240191439298911787183406185558751,    0.034319642640608095038714683012872940214,
  -0.75789963770882114801220999680989894474,   -0.045619240191439298911787183406185558751,    -0.75789963770882114801220999680989894474,    0.034319642640608095038714683012872940214,
  -0.43858148439091840506379282297401655176,    -0.75789963770882114801220999680989894474,   -0.045619240191439298911787183406185558751,    0.034319642640608095038714683012872940214,
  -0.75789963770882114801220999680989894474,    -0.43858148439091840506379282297401655176,    -0.75789963770882114801220999680989894474,    0.034319642640608095038714683012872940214,
  -0.75789963770882114801220999680989894474,    -0.75789963770882114801220999680989894474,    -0.43858148439091840506379282297401655176,    0.034319642640608095038714683012872940214,
  -0.75789963770882114801220999680989894474,   -0.045619240191439298911787183406185558751,    -0.43858148439091840506379282297401655176,    0.034319642640608095038714683012872940214,
  -0.43858148439091840506379282297401655176,    -0.75789963770882114801220999680989894474,    -0.75789963770882114801220999680989894474,    0.034319642640608095038714683012872940214,
 -0.045619240191439298911787183406185558751,    -0.43858148439091840506379282297401655176,    -0.75789963770882114801220999680989894474,    0.034319642640608095038714683012872940214,
   0.18851253896001405132314006877136059652,    -0.93444106356711465845055795933535161918,    -0.31963041182578473442202415010065735815,    0.013514495573007723718021960153355702928,
   0.18851253896001405132314006877136059652,    -0.93444106356711465845055795933535161918,    -0.93444106356711465845055795933535161918,    0.013514495573007723718021960153355702928,
  -0.93444106356711465845055795933535161918,    -0.93444106356711465845055795933535161918,     0.18851253896001405132314006877136059652,    0.013514495573007723718021960153355702928,
  -0.31963041182578473442202415010065735815,     0.18851253896001405132314006877136059652,    -0.93444106356711465845055795933535161918,    0.013514495573007723718021960153355702928,
  -0.93444106356711465845055795933535161918,    -0.31963041182578473442202415010065735815,     0.18851253896001405132314006877136059652,    0.013514495573007723718021960153355702928,
  -0.93444106356711465845055795933535161918,     0.18851253896001405132314006877136059652,    -0.93444106356711465845055795933535161918,    0.013514495573007723718021960153355702928,
  -0.31963041182578473442202415010065735815,    -0.93444106356711465845055795933535161918,     0.18851253896001405132314006877136059652,    0.013514495573007723718021960153355702928,
  -0.93444106356711465845055795933535161918,    -0.31963041182578473442202415010065735815,    -0.93444106356711465845055795933535161918,    0.013514495573007723718021960153355702928,
  -0.93444106356711465845055795933535161918,    -0.93444106356711465845055795933535161918,    -0.31963041182578473442202415010065735815,    0.013514495573007723718021960153355702928,
  -0.93444106356711465845055795933535161918,     0.18851253896001405132314006877136059652,    -0.31963041182578473442202415010065735815,    0.013514495573007723718021960153355702928,
  -0.31963041182578473442202415010065735815,    -0.93444106356711465845055795933535161918,    -0.93444106356711465845055795933535161918,    0.013514495573007723718021960153355702928,
   0.18851253896001405132314006877136059652,    -0.31963041182578473442202415010065735815,    -0.93444106356711465845055795933535161918,    0.013514495573007723718021960153355702928,
   0.60235456931668878246228336815716196359,    -0.93502943687035390432897012004314760659,    -0.73229569557598097380434312807086675041,   0.0087681963693812055566076536006009347954,
   0.60235456931668878246228336815716196359,    -0.93502943687035390432897012004314760659,    -0.93502943687035390432897012004314760659,   0.0087681963693812055566076536006009347954,
  -0.93502943687035390432897012004314760659,    -0.93502943687035390432897012004314760659,     0.60235456931668878246228336815716196359,   0.0087681963693812055566076536006009347954,
  -0.73229569557598097380434312807086675041,     0.60235456931668878246228336815716196359,    -0.93502943687035390432897012004314760659,   0.0087681963693812055566076536006009347954,
  -0.93502943687035390432897012004314760659,    -0.73229569557598097380434312807086675041,     0.60235456931668878246228336815716196359,   0.0087681963693812055566076536006009347954,
  -0.93502943687035390432897012004314760659,     0.60235456931668878246228336815716196359,    -0.93502943687035390432897012004314760659,   0.0087681963693812055566076536006009347954,
  -0.73229569557598097380434312807086675041,    -0.93502943687035390432897012004314760659,     0.60235456931668878246228336815716196359,   0.0087681963693812055566076536006009347954,
  -0.93502943687035390432897012004314760659,    -0.73229569557598097380434312807086675041,    -0.93502943687035390432897012004314760659,   0.0087681963693812055566076536006009347954,
  -0.93502943687035390432897012004314760659,    -0.93502943687035390432897012004314760659,    -0.73229569557598097380434312807086675041,   0.0087681963693812055566076536006009347954,
  -0.93502943687035390432897012004314760659,     0.60235456931668878246228336815716196359,    -0.73229569557598097380434312807086675041,   0.0087681963693812055566076536006009347954,
  -0.73229569557598097380434312807086675041,    -0.93502943687035390432897012004314760659,    -0.93502943687035390432897012004314760659,   0.0087681963693812055566076536006009347954,
   0.60235456931668878246228336815716196359,    -0.73229569557598097380434312807086675041,    -0.93502943687035390432897012004314760659,   0.0087681963693812055566076536006009347954,
    0.2561436909507320213865521444358193313,    -0.65004131563212195143010154694337920556,    -0.95606105968648811852634905054906092018,    0.017209381065149320852393906999331987612,
    0.2561436909507320213865521444358193313,    -0.65004131563212195143010154694337920556,    -0.65004131563212195143010154694337920556,    0.017209381065149320852393906999331987612,
  -0.65004131563212195143010154694337920556,    -0.65004131563212195143010154694337920556,      0.2561436909507320213865521444358193313,    0.017209381065149320852393906999331987612,
  -0.95606105968648811852634905054906092018,      0.2561436909507320213865521444358193313,    -0.65004131563212195143010154694337920556,    0.017209381065149320852393906999331987612,
  -0.65004131563212195143010154694337920556,    -0.95606105968648811852634905054906092018,      0.2561436909507320213865521444358193313,    0.017209381065149320852393906999331987612,
  -0.65004131563212195143010154694337920556,      0.2561436909507320213865521444358193313,    -0.65004131563212195143010154694337920556,    0.017209381065149320852393906999331987612,
  -0.95606105968648811852634905054906092018,    -0.65004131563212195143010154694337920556,      0.2561436909507320213865521444358193313,    0.017209381065149320852393906999331987612,
  -0.65004131563212195143010154694337920556,    -0.95606105968648811852634905054906092018,    -0.65004131563212195143010154694337920556,    0.017209381065149320852393906999331987612,
  -0.65004131563212195143010154694337920556,    -0.65004131563212195143010154694337920556,    -0.95606105968648811852634905054906092018,    0.017209381065149320852393906999331987612,
  -0.65004131563212195143010154694337920556,      0.2561436909507320213865521444358193313,    -0.95606105968648811852634905054906092018,    0.017209381065149320852393906999331987612,
  -0.95606105968648811852634905054906092018,    -0.65004131563212195143010154694337920556,    -0.65004131563212195143010154694337920556,    0.017209381065149320852393906999331987612,
    0.2561436909507320213865521444358193313,    -0.95606105968648811852634905054906092018,    -0.65004131563212195143010154694337920556,    0.017209381065149320852393906999331987612}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, simplex)] = cubDataVal(tempCoords, tempWeights);
    // order 10 and 11 orthotope
    tempNumIPs = 90;
    nIPMap[cubDataKey(tempDim, 10, orthotope)] = tempNumIPs;
    nIPMap[cubDataKey(tempDim, 11, orthotope)] = tempNumIPs;
    // nIP 90
    tempCoords.resize(tempNumIPs);
    for(int i = 0; i < tempNumIPs; i++){
      tempCoords[i].resize(tempDim);
    }
    tempWeights.resize(tempNumIPs);
    tempPvec = new std::vector<double>(
        {-0.81261433409962649639237559737974432611,                                          0,                                          0,    0.2024770736128001905853371309670196589,
0,                                          0,   0.81261433409962649639237559737974432611,    0.2024770736128001905853371309670196589,
0,   0.81261433409962649639237559737974432611,                                          0,    0.2024770736128001905853371309670196589,
0,                                          0,  -0.81261433409962649639237559737974432611,    0.2024770736128001905853371309670196589,
0.81261433409962649639237559737974432611    ,                                      0,                                          0,    0.2024770736128001905853371309670196589,
0,  -0.81261433409962649639237559737974432611,                                          0,    0.2024770736128001905853371309670196589,
0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
-0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
-0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
-0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
-0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
0.60167526419826270163441300578531749659,  -0.60167526419826270163441300578531749659,   0.60167526419826270163441300578531749659,   0.11753834795645628038993180401068212711,
0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
-0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
-0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
-0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
-0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
0.85545576101775998467509147069034657598,  -0.85545576101775998467509147069034657598,   0.85545576101775998467509147069034657598,  0.044643912078829241641001154282130043664,
0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
-0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
-0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
-0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
-0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
0.31339340451605472104577323055795129941,  -0.31339340451605472104577323055795129941,   0.31339340451605472104577323055795129941,   0.21599204525496912931346666638444131361,
-0.73466828699700801734638476986754918792,  -0.73466828699700801734638476986754918792,                                          0,   0.14519934586011569829250580079425982305,
0.73466828699700801734638476986754918792,                                          0,  -0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
0,   0.73466828699700801734638476986754918792,  -0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
0.73466828699700801734638476986754918792,   0.73466828699700801734638476986754918792,                                          0,   0.14519934586011569829250580079425982305,
0.73466828699700801734638476986754918792,                                          0,   0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
0,  -0.73466828699700801734638476986754918792,   0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
0,  -0.73466828699700801734638476986754918792,  -0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
-0.73466828699700801734638476986754918792,                                          0,   0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
-0.73466828699700801734638476986754918792,   0.73466828699700801734638476986754918792,                                          0,   0.14519934586011569829250580079425982305,
0.73466828699700801734638476986754918792,  -0.73466828699700801734638476986754918792,                                          0,   0.14519934586011569829250580079425982305,
0,   0.73466828699700801734638476986754918792,   0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
-0.73466828699700801734638476986754918792,                                          0,  -0.73466828699700801734638476986754918792,   0.14519934586011569829250580079425982305,
0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
-0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
-0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
-0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.45079993511450943037788434573026952398,   0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
-0.45079993511450943037788434573026952398,   0.96509966551271026293028182312534456821,   0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,  -0.96509966551271026293028182312534456821,  0.061441994097835335202750044633046200824,
-0.96509966551271026293028182312534456821,  -0.45079993511450943037788434573026952398,  -0.45079993511450943037788434573026952398,  0.061441994097835335202750044633046200824,
0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
-0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
-0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
-0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.94124485721060326391115015763113464139,   0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
-0.94124485721060326391115015763113464139,   0.35390281459663013491031287081289167626,   0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131,
0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,  -0.35390281459663013491031287081289167626,  0.022614296138821884223196230668984478131,
-0.35390281459663013491031287081289167626,  -0.94124485721060326391115015763113464139,  -0.94124485721060326391115015763113464139,  0.022614296138821884223196230668984478131}
);
    index = 0;
    for(int j = 0; j < tempNumIPs; j++){
      for(int i = 0; i < tempDim; i++){
        tempCoords[j][i] = (*tempPvec)[index];
        index += 1;
      }
      tempWeights[j] = (*tempPvec)[index];
      index += 1;
    }
    delete tempPvec;
    rulesDatabase[cubDataKey(tempDim, tempNumIPs, orthotope)] = cubDataVal(tempCoords, tempWeights);
  };//initializeDatabase

  bool Cubature::rulesExist = 0;
  int Cubature::maxDim;
  std::map<int, int> Cubature::maxOrderMap;
  std::map<std::tuple<int, int, elementGeometry>, int> Cubature::nIPMap;
  std::map<std::tuple<int, int, elementGeometry>,
      std::tuple< std::vector< std::vector<double> >, std::vector<double> > > Cubature::rulesDatabase;
}

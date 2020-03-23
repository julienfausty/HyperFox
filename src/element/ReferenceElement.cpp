#include "ReferenceElement.h"

namespace hfox{

  ReferenceElement::ReferenceElement(int dim, int ord, std::string geom){
    if(!databaseExists){
      initializeDatabase();
      databaseExists = 1;
    }
    setGeometry(geom);
    setDim(dim);
    setOrder(ord);
    determineNumNodes();
    determineNodes();
    determineCubature();
    determineNumFaces();
    determineFaceNodes();
    setFaceElement();
    computeInverseVandermonde();
    computeIPShapeFunctions();
    computeIPDerivShapeFunctions();
    computeDerivShapeFunctions();
  };//Constructor

  void ReferenceElement::setGeometry(std::string geom){
    if(geom == "simplex"){
      geometry = simplex;
    } else if((geom == "orthotope") or (geom == "quad") or (geom == "hex")){
      geometry = orthotope;
    } else{
      ErrorHandle eh("ReferenceElement", "setGeometry", "Element type " + geom + " is not yet supported.");
      throw(eh);
    }
  };//setGeometry

  void ReferenceElement::setDim(int dim){
    if(dim > maxDim){
      ErrorHandle eh("ReferenceElement", "setDim", "The spatial dimension " +
          std::to_string(dim) + " is not yet supported.");
      throw(eh);
    }else if(dim < 0){
      ErrorHandle eh("ReferenceElement", "setDim", "The spatial dimension " + 
          std::to_string(dim) + " is negative.");
      throw(eh);
    } else{
      dimension = dim;
    }
  };//setDim

  void ReferenceElement::setOrder(int ord){
    std::tuple<int, elementGeometry> intElem(dimension, geometry);
    int maxOrder = maxOrderMap[intElem];
    if(ord > maxOrder){
      ErrorHandle eh("ReferenceElement", "setOrder", "The interpolation order " + 
          std::to_string(ord) + " is not yet supported for dimension " + std::to_string(dimension) + ".");
      throw(eh);
    }else if(ord < 0){
      ErrorHandle eh("ReferenceElement", "setOrder", "The interpolation order " 
          + std::to_string(ord) + " is negative.");
      throw(eh);
    } else{
      order = ord;
    }
  };//setOrder

  void ReferenceElement::determineNumNodes(){
    switch(geometry){
      case simplex:
        {
          nNodes = 0;
          double factor = boost::math::factorial<double>(dimension-1);
          for(int k = 0; k < (order+1); k++){
            nNodes += 
              boost::math::factorial<double>(dimension + k - 1)/(boost::math::factorial<double>(k)*factor);
          }
          break;
        }
      case orthotope:
        nNodes = std::pow((order+1), dimension);
        break;
      default:
        ErrorHandle eh("ReferenceElement", "determineNumNodes", 
            "Something is definitely off here, if you're developping a new reference element,"
            " you forgot to write how to calculate the number of nodes.");
        throw(eh);
    }
  };//determineNumNodes

  void ReferenceElement::determineNodes(){
    nodes = nodeDatabase[std::tuple<int, int, elementGeometry>(dimension, nNodes, geometry)];
  };//determineNodes

  // Watch out here, because the finite element method needs to integrate multiplications of basis functions,
  // the cubature rules need to be higher order than just the order of the element
  void ReferenceElement::determineCubature(){
    switch(geometry){
      case simplex:
        cubatureRule = new Cubature(dimension, 2*order, geometry);
        break;
      case orthotope:
        cubatureRule = new Cubature(dimension, 4*order, geometry);
        break;
      default:
        ErrorHandle eh("ReferenceElement", "determineCubature", 
            "Something is definitely off here, if you're developping a new reference element,"
            " you forgot to write the cubature generation part.");
        throw(eh);
    }
  };//determineCubature


  void ReferenceElement::determineNumFaces(){
    switch(geometry){
      case simplex:
        nFaces = dimension + 1;
      case orthotope:
        nFaces = dimension * 2;
      default:
        ErrorHandle eh("ReferenceElement", "determineNumFaces", 
            "Something is definitely off here, if you're developping a new reference element,"
            " you forgot to write the number of faces part.");
        throw(eh);
    }
  };//determineNumFaces

  void ReferenceElement::determineFaceNodes(){
  };//determineFaceNodes

  void ReferenceElement::setFaceElement(){
    std::map<elementGeometry, std::string> e2stringMap;
    e2stringMap[simplex] = "simplex";
    e2stringMap[orthotope] = "orthotope";
    if(dimension != 0){
      faceElement = new ReferenceElement(dimension - 1, order, e2stringMap[geometry]);
    } else{
      faceElement = NULL;
    }
  };//setFaceElement

  void ReferenceElement::computeInverseVandermonde(){
    EMatrix V = EMatrix::Constant(nNodes, nNodes, 1.0);
    std::vector< std::vector<int> > nodeToModeMap;
    std::vector<int> orders(order+1);
    std::iota(std::begin(orders), std::end(orders), 0);
    nodeToModeMap = generateCombinationsWithRepetition(orders, dimension);
    std::vector<int> mode;
    double coord;
    using boost::math::jacobi;
    switch(geometry){
      case orthotope:
        {
          //mode loop
          for(int i = 0; i < nNodes; i++){
            mode = nodeToModeMap[i];
            //node loop
            for(int j = 0; j < nNodes; j++){
              //dim loop
              for(int k = 0; k < dimension; k++){
                V(j, i) = jacobi(mode[k], 0.0, 0.0, nodes[j][k]);        
              }
            }
          }
          break;
        }
      case simplex:
        {
          std::vector< std::vector<int> > realNodeToModeMap;
          for(int i = 0; i < nodeToModeMap.size(); i++){
            int sum = std::accumulate(nodeToModeMap[i].begin(), nodeToModeMap[i].end(), 0);
            if(sum <= order){
              realNodeToModeMap.push_back(nodeToModeMap[i]);
            }
          }
          std::vector<double> mappedCoords(dimension);
          //node loop
          for(int i = 0; i < nNodes; i++){
            switch(dimension){
              case 1:
                {
                  mappedCoords = nodes[i];
                  break;
                }
              case 2:
                {
                  mappedCoords[0] = 2.0*(1.0+nodes[i][0])/(1.0-nodes[i][1]) - 1.0;
                  mappedCoords[1] = nodes[i][1];
                  break;
                }
              case 3:
                {
                  mappedCoords[0] = -2.0*(1.0+nodes[i][0])/(nodes[i][1] + nodes[i][2])-1.0;
                  mappedCoords[1] = 2.0*(1.0 + nodes[i][1])/(1.0-nodes[i][2]) - 1.0;
                  mappedCoords[2] = nodes[i][2];
                  break;
                }
              }
            //mode loop
            for(int j = 0; j < nNodes; j++){
              mode = realNodeToModeMap[j];
              //dim loop
              for(int k = 0; k < dimension; k++){
                switch(k){
                  case 0:
                    {
                      V(i, j) *= jacobi(mode[k], 0.0, 0.0, mappedCoords[k]);
                      break;
                    }
                  case 1:
                    {
                      V(i, j) *= std::sqrt(2.0)*jacobi(mode[k], 2.0*mode[k-1] + 1.0, 0.0, mappedCoords[k])*std::pow(1-mappedCoords[k], mode[k-1]);
                      break;
                    }
                  case 2:
                    {
                      V(i, j) *= 2.0*jacobi(mode[k], 2.0*(mode[k-1] + mode[k-2] + 1.0), 0.0, mappedCoords[k])*std::pow(1-mappedCoords[k], mode[k-1]+mode[k-2]);
                    break;
                    }
                }
              }
            }
            }
          break;
        }
    }
    invV = V.inverse();
  };//computeInverseVandermonde

  void ReferenceElement::computeIPShapeFunctions(){
    int nIPs = cubatureRule->getNumIPs();
    const std::vector< std::vector<double> > * IPs = cubatureRule->getIPCoords();
    ipShapeFunctions.resize(nIPs);
    for(int i = 0; i < nIPs; i++){
      ipShapeFunctions[i] = interpolate((*IPs)[i]);
    }
  };//computeIPShapeFunctions

  void ReferenceElement::computeIPDerivShapeFunctions(){
    int nIPs = cubatureRule->getNumIPs();
    const std::vector< std::vector<double> > * IPs = cubatureRule->getIPCoords();
    ipDerivShapeFunctions.resize(nIPs);
    for(int i = 0; i < nIPs; i++){
      ipDerivShapeFunctions[i] = interpolateDeriv((*IPs)[i], 1);
    }
  };//computeIPDerivShapeFunctions

  void ReferenceElement::computeDerivShapeFunctions(){
    derivShapeFunctions.resize(nNodes);
    for(int i = 0; i < nNodes; i++){
      derivShapeFunctions[i] = interpolateDeriv(nodes[i], 1);
    }
  };//computeDerivShapeFunctions

  int ReferenceElement::getDimension() const{
    return dimension;
  };//getDimension

  int ReferenceElement::getOrder() const{
    return order;
  };//getOrder

  elementGeometry ReferenceElement::getGeometry() const{
    return geometry;
  };//getGeometry

  int ReferenceElement::getNumNodes() const{
    return nNodes;
  };//getNumNodes

  int ReferenceElement::getNumIPs() const{
    return cubatureRule->getNumIPs();
  };//getNumIPs

  int ReferenceElement::getNumFaces() const{
    return nFaces;
  };//getNumFaces

  const std::vector< std::vector<double> > * ReferenceElement::getNodes() const{
    return &nodes;
  };//getNodes

  const std::vector< std::vector<int> > * ReferenceElement::getFaceNodes() const{
    return &faceNodes;
  };//getFaceNodes

  const std::vector<int> * ReferenceElement::getInnerNodes() const{
    return &innerNodes;
  };//getInnerNodes

  const std::vector< std::vector<double> > * ReferenceElement::getIPCoords() const{
    return cubatureRule->getIPCoords();
  };//getIPCoords

  const std::vector<double> * ReferenceElement::getIPWeights() const{
    return cubatureRule->getIPWeights();
  };//getIPWeights

  const std::vector< std::vector<double> > * ReferenceElement::getIPShapeFunctions() const{
    return &ipShapeFunctions;
  };//getIPShapeFunctions

  const std::vector< std::vector< std::vector<double> > > * ReferenceElement::getIPDerivShapeFunctions() const{
    return &ipDerivShapeFunctions;
  };//getIPDerivShapeFunctions

  const std::vector< std::vector< std::vector<double> > > * ReferenceElement::getDerivShapeFunctions() const{
    return &derivShapeFunctions;
  };//getDerivShapeFunctions

  const ReferenceElement * ReferenceElement::getFaceElement() const{
    return faceElement;
  };//getFaceElement

  std::vector<double> ReferenceElement::interpolate(const std::vector<double> & point) const{
  };//interpolate

  ReferenceElement::~ReferenceElement(){
    delete cubatureRule;
    if(dimension != 0){
      delete faceElement;
    }
  };//Destructor

  std::vector< std::vector<int> > ReferenceElement::generateCombinationsWithRepetition(std::vector<int> & set, int size) const{
    std::vector< std::vector<int> > combinations;
    if(size == 1){
      combinations.resize(set.size());
      for(int i = 0; i < set.size(); i++){
        combinations[i] = std::vector<int>(1, set[i]);
      }
      return combinations;
    } else if(size > 1){
      std::vector< std::vector<int> > iCombinations(set.size());
      std::vector<int> newcombination(size);
      for(int i = 0; i < set.size(); i++){
        iCombinations = generateCombinationsWithRepetition(set, size-1);
        newcombination[0] = set[i];
        for(int j = 0; j < iCombinations.size(); j++){
          for(int k = 0; k < iCombinations[j].size(); k++){
            newcombination[k+1] = iCombinations[j][k];
          }
        }
        combinations.push_back(newcombination);
      }
      return combinations;
    } else{
      ErrorHandle eh("ReferenceElement", "generateCombinationsWithRepetition", "the size passed to "
          "generateCombinationWithRepetion is lower than 1.");
      throw(eh);
    }
  }

  //static variables
  bool databaseExists = 0;
  int maxDim = 3;
  std::map<std::tuple<int, elementGeometry>, int> maxOrderMap;
  std::map< std::tuple<int, int, elementGeometry>, std::vector< std::vector<double> > > nodeDatabase;

}//hfox

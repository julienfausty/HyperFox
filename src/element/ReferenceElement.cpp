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
    if(dimension > 0){
      determineNodeToModeMap();
      computeInverseVandermonde();
      computeIPShapeFunctions();
      computeIPDerivShapeFunctions();
      computeDerivShapeFunctions();
    }
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
    if(dimension > 0){
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
    }else{
      nNodes = 0;
    }
  };//determineNumNodes

  void ReferenceElement::determineNodes(){
    nodes = nodeDatabase[std::tuple<int, int, elementGeometry>(dimension, order, geometry)];
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
        break;
      case orthotope:
        nFaces = dimension * 2;
        break;
      default:
        ErrorHandle eh("ReferenceElement", "determineNumFaces", 
            "Something is definitely off here, if you're developping a new reference element,"
            " you forgot to write the number of faces part.");
        throw(eh);
    }
  };//determineNumFaces

  void ReferenceElement::determineFaceNodes(){
  };//determineFaceNodes

  void ReferenceElement::determineNodeToModeMap(){
    std::vector<int> orders(order+1);
    std::iota(std::begin(orders), std::end(orders), 0);
    nodeToModeMap = generateCombinationsWithRepetition(orders, dimension);
    if(geometry == simplex){
      std::vector< std::vector<int> > realNodeToModeMap;
      for(int i = 0; i < nodeToModeMap.size(); i++){
        int sum = std::accumulate(nodeToModeMap[i].begin(), nodeToModeMap[i].end(), 0);
        if(sum <= order){
          realNodeToModeMap.push_back(nodeToModeMap[i]);
        }
      }
      nodeToModeMap = realNodeToModeMap;
    }
  };//determineNodeToModeMap

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


  std::vector<double> ReferenceElement::mapCoordsSimplex(const std::vector<double> & coords) const{
    std::vector<double> mappedCoords(dimension);
    switch(dimension){
      case 1:
        {
          mappedCoords = coords;
          break;
        }
      case 2:
        {
          if(coords[1] != 1.0){
            mappedCoords[0] = 2.0*(1.0+coords[0])/(1.0-coords[1]) - 1.0;
          } else{
            mappedCoords[0] = -1.0;
          }
          mappedCoords[1] = coords[1];
          break;
        }
      case 3:
        {
          if((coords[1] + coords[2]) != 0.0){
            mappedCoords[0] = -2.0*(1.0+coords[0])/(coords[1] + coords[2])-1.0;
          } else {
            mappedCoords[0] = -1.0;
          }
          if(coords[2] != 1.0){
            mappedCoords[1] = 2.0*(1.0 + coords[1])/(1.0-coords[2]) - 1.0;
          } else{
            mappedCoords[1] = -1.0;
          }
          mappedCoords[2] = coords[2];
          break;
        }
    }
    return mappedCoords;
  };//mapCoordsSimplex

  std::vector<double> ReferenceElement::computeModes(const std::vector<double> & point) const{
    std::vector<double> result(nNodes, 1.0);
    std::vector<int> mode;
    using boost::math::jacobi;
    switch(geometry){
      case orthotope:
        {
          //mode loop
          for(int i = 0; i < nNodes; i++){
            mode = nodeToModeMap[i];
            //dim loop
            for(int k = 0; k < dimension; k++){
              result[i] *= jacobi(mode[k], 0.0, 0.0, point[k]);        
            }
          }
          break;
        }
      case simplex:
        {
          std::vector<double> mappedCoords = mapCoordsSimplex(point);
          //mode loop
          for(int j = 0; j < nNodes; j++){
            mode = nodeToModeMap[j];
            //dim loop
            for(int k = 0; k < dimension; k++){
              switch(k){
                case 0:
                  {
                    result[j] *= jacobi(mode[k], 0.0, 0.0, mappedCoords[k]);
                    break;
                  }
                case 1:
                  {
                    result[j] *= std::sqrt(2.0)*jacobi(mode[k], 2.0*mode[k-1] + 1.0, 0.0, mappedCoords[k])*std::pow(1-mappedCoords[k], mode[k-1]);
                    break;
                  }
                case 2:
                  {
                    result[j] *= 2.0*jacobi(mode[k], 2.0*(mode[k-1] + mode[k-2] + 1.0), 0.0, mappedCoords[k])*std::pow(1-mappedCoords[k], mode[k-1]+mode[k-2]);
                    break;
                  }
              }
            }
          }
          break;
        }
    }
    return result;
  };//computeModes


  std::vector< std::vector<double> > ReferenceElement::computeDerivModes(const std::vector<double> & point) const{
    std::vector< std::vector<double> > result(nNodes, std::vector<double>(dimension, 1.0));
    std::vector<int> mode;
    using boost::math::jacobi;
    using boost::math::jacobi_derivative;
    switch(geometry){
      case orthotope:
        {
          //mode loop
          for(int i = 0; i < nNodes; i++){
            mode = nodeToModeMap[i];
            //dim loops
            for(int k = 0; k < dimension; k++){
              for(int l = 0; l < dimension; l++){
                if(l != k){
                  result[i][l] *= jacobi(mode[k], 0.0, 0.0, point[k]); 
                } else{
                  result[i][l] *= jacobi_derivative(mode[k], 0.0, 0.0, point[k], 1);
                }
              }
            }
          }
          break;
        }
      case simplex:
        {
          std::vector<double> mappedCoords = mapCoordsSimplex(point);
          EMatrix Jacobian(dimension, dimension);
          switch(dimension){
            case 1:
              {
                Jacobian(0,0) = 1.0;
                break;
              }
            case 2:
              {
                if(point[1] != 1.0){
                  Jacobian(0, 0) = 2.0/(1.0-point[1]); 
                  Jacobian(1, 0) = 2.0*(1.0 + point[0])*(1.0/std::pow(1.0-point[1], 2.0));
                }else{
                  Jacobian(0, 0) = 2.0;
                  Jacobian(1, 0) = 0.0;
                }
                Jacobian(0, 1) = 0.0;
                Jacobian(1,1) = 1.0;
                break;
              }
            case 3:
              {
                if((point[1] + point[2]) != 0.0){
                  Jacobian(0,0) = -2.0/(point[1]+point[2]);
                  Jacobian(1, 0) = 2.0*(1.0+point[0])*(1.0/std::pow(point[1] + point[2], 2.0));
                  Jacobian(2, 0) = 2.0*(1.0+point[0])*(1.0/std::pow(point[1] + point[2], 2.0));
                } else{
                  Jacobian(0, 0) = 2.0;
                  Jacobian(1, 0) = 0.0;
                  Jacobian(2, 0) = 0.0;
                }
                if(point[2] != 1.0){
                  Jacobian(1, 1) = 2.0/(1.0 - point[2]); 
                  Jacobian(2, 1) = 2.0*(1.0 + point[1])*(1.0/std::pow(1.0-point[2], 2.0));
                }else{
                  Jacobian(1, 1) = 2.0;
                  Jacobian(2, 1) = 0.0;
                }
                Jacobian(0, 1) = 0.0;
                Jacobian(0, 2) = 0.0;
                Jacobian(1, 2) = 0.0;
                Jacobian(2, 2) = 1.0;
                break;
              }
          }
          EVector gradMode(dimension);
          //mode loop
          for(int j = 0; j < nNodes; j++){
            mode = nodeToModeMap[j];
            gradMode = EVector::Constant(dimension, 1.0);
            //dim loop
            for(int k = 0; k < dimension; k++){
              for(int l = 0; l < dimension; l++){
                switch(l){
                  case 0:
                    {
                      if(l!=k){
                        gradMode[k] *= jacobi(mode[l], 0.0, 0.0, mappedCoords[l]);
                      }else{
                        gradMode[k] *= jacobi_derivative(mode[l], 0.0, 0.0, mappedCoords[l], 1);
                      }
                      break;
                    }
                  case 1:
                    {
                      if(l != k){
                        gradMode[k] *= std::sqrt(2.0)*jacobi(mode[l], 2.0*mode[l-1] + 1.0, 0.0, mappedCoords[l])*std::pow(1-mappedCoords[l], mode[l-1]);
                      } else{
                        int power = mode[l-1] - 1;
                        if (power < 0){ power = 0;}
                        gradMode[k] *= std::sqrt(2.0)*jacobi_derivative(mode[l], 2.0*mode[l-1] + 1.0, 0.0, mappedCoords[l], 1)*std::pow(1-mappedCoords[l], mode[l-1]) +
                          std::sqrt(2.0)*jacobi(mode[l], 2.0*mode[l-1] + 1.0, 0.0, mappedCoords[l])*(-mode[l-1]*std::pow(1-mappedCoords[l], power));
                      }
                      break;
                    }
                  case 2:
                    {
                      if(k != l){
                        gradMode[k] *= 2.0*jacobi(mode[l], 2.0*(mode[l-1] + mode[l-2] + 1.0), 0.0, mappedCoords[l])*std::pow(1-mappedCoords[l], mode[l-1]+mode[l-2]);
                      }else{
                        int power = mode[l-1] + mode[l-2] - 1;
                        if (power < 0){ power = 0;}
                        gradMode[k] *= 2.0*jacobi_derivative(mode[l], 2.0*(mode[l-1] + mode[l-2] + 1.0), 0.0, mappedCoords[l], 1)*std::pow(1-mappedCoords[l], mode[l-1]+mode[l-2]) +
                          2.0*jacobi(mode[l], 2.0*(mode[l-1] + mode[l-2] + 1.0), 0.0, mappedCoords[l])*(-(mode[l-1] + mode[l-2])*std::pow(1-mappedCoords[l], power));
                      }
                      break;
                    }
                }
              }
            }
            gradMode = Jacobian * gradMode;
            result[j] = std::vector<double>(gradMode.data(), gradMode.data() + dimension);
          }
          break;
        }
    }
    return result;
  };//computeModes

  void ReferenceElement::computeInverseVandermonde(){
    EMatrix V(nNodes, nNodes);
    std::vector<double> oneRow(nNodes);
    for(int i = 0; i < nNodes; i++){
      oneRow = computeModes(nodes[i]);
      V.row(i) = EVector::Map(oneRow.data(), nNodes);
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
      ipDerivShapeFunctions[i] = interpolateDeriv((*IPs)[i]);
    }
  };//computeIPDerivShapeFunctions

  void ReferenceElement::computeDerivShapeFunctions(){
    derivShapeFunctions.resize(nNodes);
    for(int i = 0; i < nNodes; i++){
      derivShapeFunctions[i] = interpolateDeriv(nodes[i]);
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
    std::vector<double> modesAtPoint = computeModes(point);
    EVector valModes(nNodes);
    valModes = EVector::Map(modesAtPoint.data(), nNodes);
    EMatrix lagrangePolys(1, nNodes);
    lagrangePolys = valModes.transpose() * invV;
    std::vector<double> result(lagrangePolys.data(), lagrangePolys.data() + nNodes);
    return result;
  };//interpolate

  std::vector< std::vector<double> > ReferenceElement::interpolateDeriv(const std::vector<double> & point) const{
    std::vector< std::vector<double> > modesAtPoint = computeDerivModes(point);
    std::vector< std::vector<double> > result(nNodes, std::vector<double>(dimension));
    EVector valModes(nNodes);
    EMatrix lagrangePolys(1, nNodes);
    for(int k = 0; k < dimension; k++){
      //mode loop
      for(int i = 0; i < nNodes; i++){
        valModes[i] = modesAtPoint[i][k];
      }
      lagrangePolys = valModes.transpose() * invV;
      for(int i = 0; i < nNodes; i++){
        result[i][k] = lagrangePolys(0, i);
      }
    }
    return result;
  };//interpolateDeriv

  ReferenceElement::~ReferenceElement(){
    delete cubatureRule;
    if(dimension != 0){
      delete faceElement;
    }
  };//Destructor

  std::vector< std::vector<int> > ReferenceElement::generateCombinationsWithRepetition(std::vector<int> & set, int size){
    std::vector< std::vector<int> > combinations;
    if(size == 1){
      combinations.resize(set.size());
      for(int i = 0; i < set.size(); i++){
        combinations[i] = std::vector<int>(1, set[i]);
      }
      return combinations;
    } else if(size > 1){
      std::vector< std::vector<int> > iCombinations = generateCombinationsWithRepetition(set, size-1);
      std::vector<int> newcombination(size);
      for(int i = 0; i < set.size(); i++){
        newcombination[0] = set[i];
        for(int j = 0; j < iCombinations.size(); j++){
          for(int k = 0; k < iCombinations[j].size(); k++){
            newcombination[k+1] = iCombinations[j][k];
          }
          combinations.push_back(newcombination);
        }
      }
      return combinations;
    } else{
      ErrorHandle eh("ReferenceElement", "generateCombinationsWithRepetition", "the size passed to "
          "generateCombinationWithRepetion is lower than 1.");
      throw(eh);
    }
  }

  void ReferenceElement::initializeDatabase(){
    typedef std::tuple<int, elementGeometry> maxOrderKey;
    typedef std::tuple<int, int, elementGeometry> dataKey;
    typedef std::vector< std::vector<double> > dataVal;

    maxDim = 3;
    maxOrderMap[maxOrderKey(0, simplex)] = 10;
    maxOrderMap[maxOrderKey(1, simplex)] = 10;
    maxOrderMap[maxOrderKey(2, simplex)] = 10;
    maxOrderMap[maxOrderKey(3, simplex)] = 5;
    maxOrderMap[maxOrderKey(0, orthotope)] = 5;
    maxOrderMap[maxOrderKey(1, orthotope)] = 5;
    maxOrderMap[maxOrderKey(2, orthotope)] = 5;
    maxOrderMap[maxOrderKey(3, orthotope)] = 2;

    int maxOrder0 = std::max(maxOrderMap[maxOrderKey(0, simplex)], maxOrderMap[maxOrderKey(0, orthotope)]);
    for(int i = 0; i < maxOrder0; i++){
      nodeDatabase[dataKey(0, i, simplex)] = dataVal();
      nodeDatabase[dataKey(0, i, orthotope)] = dataVal();
    }
    int tempDim, tempOrder, tempNNodes, index;
    dataVal tempVal;
    std::vector<double> * tempVec;
    //orthotope and simplex dimension 1
    tempDim = 1;
    // order 0
    tempOrder = 0;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { 0.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 1
    tempOrder = 1;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 2
    tempOrder = 2;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, 0.0, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 3
    tempOrder = 3;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.447213595499957939282, 0.447213595499957939282, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 4
    tempOrder = 4;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.6546536707079771437983, 0.0, 0.6546536707079771437983, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 5
    tempOrder = 5;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.765055323929464692851,
        -0.2852315164806450963142, 0.2852315164806450963142,
        0.765055323929464692851, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 6
    tempOrder = 6;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.830223896278566929872,
        -0.4688487934707142138038, 0.0, 0.4688487934707142138038,
        0.830223896278566929872, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 7
    tempOrder = 7;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.8717401485096066153375,
        -0.5917001814331423021445, -0.2092992179024788687687, 
        0.2092992179024788687687, 0.5917001814331423021445, 
        0.8717401485096066153375, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 8
    tempOrder = 8;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.8997579954114601573124,
        -0.6771862795107377534459, -0.3631174638261781587108,
        0.0,
        0.3631174638261781587108, 0.6771862795107377534459, 
        0.8997579954114601573124, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 9
    tempOrder = 9;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.9195339081664588138289,
        -0.7387738651055050750031, -0.4779249498104444956612,
        -0.1652789576663870246262, 0.1652789576663870246262,
        0.4779249498104444956612, 0.7387738651055050750031, 
        0.9195339081664588138289, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;
    // order 10
    tempOrder = 10;
    tempNNodes = tempOrder + 1;
    tempVal.resize(tempNNodes);
    for(int i = 0; i < tempNNodes; i++){
      tempVal[i].resize(tempDim);
    }
    tempVec = new std::vector<double>(
        { -1.0, -0.9340014304080591343323,
        -0.7844834736631444186224, -0.565235326996205006471,
        -0.2957581355869393914319, 0.0, 0.2957581355869393914319,
        0.565235326996205006471, 0.7844834736631444186224, 
        0.9340014304080591343323, 1.0 }
        );
    index = 0;
    for(int i = 0; i < tempNNodes; i++){
      for(int j = 0; j < tempDim; j++){
        tempVal[i][j] = (*tempVec)[index];
        index += 1;
      }
    }
    nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;   
    nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
    delete tempVec;

    std::vector<double> lobattoPts;
    std::vector<std::vector<int> > combinations;
    //orthotope dim k
    for(int i = 2; i < (maxDim+1); i++){
      tempDim = i;
      for(int j = 0; j < (maxOrderMap[maxOrderKey(i, orthotope)] + 1); j++){
        tempOrder = j;
        tempNNodes = std::pow((tempOrder + 1), tempDim);
        tempVal.resize(tempNNodes);
        for(int k = 0; k < tempNNodes; k++){
          tempVal[k].resize(tempDim);
        }
        std::vector<int> indexes(tempOrder + 1);
        std::iota(indexes.begin(), indexes.end(), 0);
        combinations = generateCombinationsWithRepetition(indexes, tempDim);
        lobattoPts.resize(tempOrder + 1);
        for(int k = 0; k < tempOrder + 1; k++){
          lobattoPts[k] = (nodeDatabase[dataKey(1, tempOrder, orthotope)])[k][0];
        }
        for(int k = 0; k < combinations.size(); k++){
          for(int l = 0; l < tempDim; l++){
            tempVal[k][l] = lobattoPts[combinations[k][l]];
          }
        }
        nodeDatabase[dataKey(tempDim, tempOrder, orthotope)] = tempVal;
      }
    }

    //simplex dim k (Lobatto grid)
    for(int i = 2; i < (maxDim+1); i++){
      tempDim = i;
      for(int j = 0; j < (maxOrderMap[maxOrderKey(i, simplex)] + 1); j++){
        tempOrder = j;
        lobattoPts.resize(tempOrder + 1);
        for(int k = 0; k < tempOrder + 1; k++){
          lobattoPts[k] = (nodeDatabase[dataKey(1, tempOrder, simplex)])[k][0];
        }
        std::vector<int> indexes(tempOrder + 1);
        std::iota(indexes.begin(), indexes.end(), 0);
        combinations = generateCombinationsWithRepetition(indexes, tempDim);
        std::vector< std::vector<int> > reduced;
        int sum = 0.0;
        for(int k = 0; k < combinations.size(); k++){
          sum = std::accumulate(combinations[k].begin(), combinations[k].end(), 0);
          if(sum <= tempOrder){
            reduced.push_back(combinations[k]);
          }
        }
        combinations = reduced;
        tempNNodes = combinations.size();
        std::map<std::vector<int>, double> coefficients;
        for(int k = 0; k < tempNNodes; k++){
          coefficients[combinations[k]] = 2.0 + tempDim*lobattoPts[combinations[k][0]];
          for(int l = 1; l < tempDim; l++){
            coefficients[combinations[k]] -= lobattoPts[combinations[k][l]];
          }
          coefficients[combinations[k]] -= lobattoPts[tempOrder - std::accumulate(combinations[k].begin(), combinations[k].end(), 0)];
          coefficients[combinations[k]] *= 1.0/(tempDim + 1.0);
          coefficients[combinations[k]] -= 1.0;
        }
        tempVal.resize(tempNNodes);
        for(int k = 0; k < tempNNodes; k++){
          tempVal[k].resize(tempDim);
          switch(tempDim){
            case 2:
              {
                tempVal[k][0] = coefficients[combinations[k]];
                tempVal[k][1] = coefficients[std::vector<int>({combinations[k][1], combinations[k][0]})];
                break;
              }
            case 3:
              {
                tempVal[k][0] = coefficients[combinations[k]];
                tempVal[k][1] = coefficients[std::vector<int>({combinations[k][1], combinations[k][0], combinations[k][2]})];
                tempVal[k][2] = coefficients[std::vector<int>({combinations[k][2], combinations[k][1], combinations[k][0]})];
                break;
              }
          }
        }
        nodeDatabase[dataKey(tempDim, tempOrder, simplex)] = tempVal;
      }
    }
  };//initializeDatabase

  //static variables
  bool ReferenceElement::databaseExists = 0;
  int ReferenceElement::maxDim;
  std::map<std::tuple<int, elementGeometry>, int> ReferenceElement::maxOrderMap;
  std::map< std::tuple<int, int, elementGeometry>, std::vector< std::vector<double> > > ReferenceElement::nodeDatabase;

}//hfox

#include "ErrorHandle.h"

namespace hfox{

const char * ErrorHandle::what() const noexcept{
  std::string delim(" : ");
  std::string full_msg = userclass + delim + userfunction + delim + err_msg;
  return full_msg.c_str();
};

void ErrorHandle::checkPointList(const std::vector< std::vector<double> > & 
    points) {
  int dim = (points[0]).size();
  int currentdim;
  for(auto vit = (points.begin()+1); vit != points.end(); vit++){
    currentdim = (*vit).size();
    if(dim != currentdim){
      err_msg.assign("point list dimensions not respected! (" + 
        std::to_string(dim) + " != "+ std::to_string(currentdim) + ")\n");
      throw (*this);
      break;
    }
  }
};

void ErrorHandle::checkIndexList(int max_index, 
    const std::vector< std::vector<int> > & points) {
  for(auto vit = (points.begin()+1); vit != points.end(); vit++){
    for(auto subvit = (vit->begin()+1); subvit != vit->end(); subvit++){
      if((*subvit)  > max_index){
        err_msg.assign("index list range not respected! (" + 
            std::to_string(*subvit) + " > "+ std::to_string(max_index) + 
            ")\n");
        throw (*this);
        break;
      }
    }
  }
};

} //hfox

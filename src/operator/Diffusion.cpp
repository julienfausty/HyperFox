#include "Diffusion.h"

namespace hfox{

void Diffusion::assemble(const std::vector< double > & detJacobians, 
    const std::vector< EMatrix > & invJacobians){
  if(allocated){
  } else {
    ErrorHandle("Diffusion", "assemble", "cannot assemble before allocating.");
  }
};//assemble

}

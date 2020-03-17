#include "Io.h"

namespace hfox{

  std::string Io::getExtension(std::string filename){

      // Create a Path object from given string
    boost::filesystem::path pathObj(filename);
      // Check if file name in the path object has extension
      if (pathObj.has_extension()) {
      // Fetch the extension from path object and return
        return pathObj.extension().string();
      }
      // In case of no extension return empty string
      return "";
  }//getExtension

}//hfox

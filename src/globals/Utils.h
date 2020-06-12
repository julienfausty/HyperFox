#ifndef UTILS_H
#define UTILS_H

namespace hfox{

class Utils{

  public:
    template<class IndexContainer, class Container>
      static void slice(const IndexContainer & indexes, const Container * list2Slice, Container * slice){
        slice->resize(indexes.size());
        for(int i = 0; i < indexes.size(); i++){
          slice->at(i) = list2Slice->at(indexes[i]);
        }
      };

};//Utils

}//hfox

#endif//UTILS_H

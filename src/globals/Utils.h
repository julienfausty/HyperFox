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
    /*!
     * \brief helper function for multiplying indexes
     * @param dim number of objects per idex
     * @param idexes the list of indexes
     * @param multiIndexes the list of multiplied indexes
     */
    static void multiplyIndexes(int dim, const std::vector<int> * indexes, std::vector<int> * multiIndexes){
      multiIndexes->resize(dim*(indexes->size()));
      for(int i = 0; i < indexes->size(); i++){
        for(int j = 0; j < dim; j++){
          multiIndexes->at(i*dim + j) = indexes->at(i)*dim + j;
        }
      }
    };

};//Utils

}//hfox

#endif//UTILS_H

#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>
#include <numeric>

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

    /*!
     * \brief helper function to Allgather a set of local index lists into a global one
     * @param locList the local list
     * @param globList the global list
     * @param comm the MPI communicator
     */
    static void allGather(std::vector<int> * locList, std::vector<int> * globList, MPI_Comm comm){
      int nParts;
      MPI_Comm_size(comm, &nParts);
      int locSize = locList->size();
      std::vector<int> globSizes(nParts);
      std::vector<int> offsets(nParts);
      MPI_Allgather(&locSize, 1, MPI_INT, globSizes.data(), 1, MPI_INT, comm);
      for(int i = 0; i < nParts; i++){
        offsets[i] = std::accumulate(globSizes.begin(), globSizes.begin() + i, 0);
      }
      globList->resize(std::accumulate(globSizes.begin(), globSizes.end(), 0));
      MPI_Allgatherv(locList->data(), locList->size(), MPI_INT, globList->data(), globSizes.data(), offsets.data(), MPI_INT, comm);
    };

};//Utils

}//hfox

#endif//UTILS_H

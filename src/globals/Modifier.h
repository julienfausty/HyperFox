#ifndef MODIFIER_H
#define MODIFIER_H

#include <algorithm>

namespace hfox{

enum ModType{
  APPEND,
  REMOVE
};

/*!
 * \brief a template class that can modify stl containers in certain ways
 */

template <class Container>
class Modifier{
  public:
    /*!
     * \brief empty constructor
     */
    Modifier(){};
    /*!
     * \brief empty destructor
     */
    ~Modifier(){};
    /*!
     * \brief sets the modification container
     * @param values a container of modification values
     */
    void setValues(Container * values){
      vals = values;
    };
    /*!
     * \brief sets the modification type
     * @param type the modification type
     */
    void setType(ModType type){
      myType = type;
    };
    /*!
     * \brief modifies the container
     * @param c container to modify
     */
    template <class ModContainer>
    void modify(ModContainer * c) const{
      switch(myType){
        case APPEND:{
                      append(c);
                      break;
                    }
        case REMOVE:{
                      remove(c);
                      break;
                    }
      }
    };
    /*!
     * \brief append vals to c
     * @param c the container to append to
     */
    template <class ModContainer>
    void append(ModContainer * c) const{
      c->insert(c->end(), vals->begin(), vals->end());
    };
    /*!
     * \brief remove elements at indexes of vals from c
     * @param c the container to remove from
     */
    template <class ModContainer>
    void remove(ModContainer * c) const{
      std::sort(vals->begin(), vals->end());
      typename ModContainer::iterator readit, writeit;
      writeit = c->begin() + *(vals->begin());
      typename Container::iterator delIt = vals->begin() + 1;
      for(readit = writeit + 1; readit != c->end(); readit++){
        while((std::distance(c->begin(), readit) == *(delIt)) and (delIt != vals->end())){
          readit += 1;
          delIt += 1;
        }
        if(readit != c->end()){
          *writeit = *readit;
          writeit += 1;
        } else {
          break;
        }
      }
      c->erase(writeit, c->end());
      //ModContainer buffer(c->size() - vals->size());
      //int pos = 0;
      //for(int i = 0; i < c->size(); i++){
        //if(std::find(vals->begin(), vals->end(), i) == vals->end()){
          //buffer[pos] = c->at(i);
          //pos += 1;
        //}
      //}
      //c->resize(buffer.size());
      //std::copy(buffer.begin(), buffer.end(), c->begin());
    };
  protected:
    /*!
     * \brief the container one uses to modify
     */
    Container * vals = NULL;
    /*!
     * \brief the modification type
     */
    ModType myType;

};//Modifier

}//hfox

#endif//MODIFIER_H

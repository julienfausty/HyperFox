#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <chrono>
#include <iostream>
#include <sys/ioctl.h>
#include <unistd.h>

/*!
 * \brief A (hopefully) easy to use progress bar
 */

class ProgressBar{
  public:
    /*!
     * \brief empty constructor
     */
    ProgressBar(){};
    /*!
     * \brief function for updating the bar
     */
    void update();
    /*!
     * \brief set the iteration index
     *
     * @param pi pointer to index
     */
    void setIterIndex(int * pi){index = pi;};
    /*!
     * \brief set the total number of iterations
     */
    void setNumIterations(int numIters){nIters = numIters;};
  protected:
    /*!
     * \brief method to update display
     */
    void updateDisplay();
    /*!
     * \brief pointer to the running index
     */
    int * index;
    /*!
     * \brief number of total iterations
     */
    int nIters;
    /*!
     * \brief update frequency (i.e. update every freq seconds)
     */
    std::chrono::duration<double> freq = std::chrono::seconds(1);
    /*!
     * \brief the time of the last update
     */
    std::chrono::time_point<std::chrono::steady_clock> lastUpdate;
    /*!
     * \brief the char to use for the progress
     */
    std::string progressChar = "|";
};//ProgressBar

#endif//PROGRESSBAR_H



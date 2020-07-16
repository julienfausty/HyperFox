#include "ProgressBar.h"


void ProgressBar::setNumIterations(int numIters){
  locIters = numIters;
  MPI_Reduce(&numIters, &nIters, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
};//setNumIterations


void ProgressBar::update(){
  std::chrono::time_point<std::chrono::steady_clock> thisUpdate = std::chrono::steady_clock::now();
  if((thisUpdate - lastUpdate > freq) or (*index == (locIters-1))){
    lastUpdate = thisUpdate;
    totalIndex = *index;
    //MPI_Reduce(index, &totalIndex, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0){
      updateDisplay();
    }
  }
};//update;


void ProgressBar::updateDisplay(){
  int cols = 80;

//#ifdef TIOCGSIZE
  //struct ttysize ts;
  //ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
  //cols = ts.ts_cols;
//#elif defined(TIOCGWINSZ)
  //struct winsize ts;
  //ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
  //cols = ts.ws_col; 
//#endif [> TIOCGSIZE <]

  cols -= 4;
  int perc = ((*index + 1)*100)/(locIters);
  int nProg = (cols*perc)/100;
  std::cout << "\r" << std::string(nProg, progressChar.at(0)) << std::string(cols-nProg, ' ') << perc << "%" << std::flush;
  if(perc == 100){
    std::cout << std::endl;
  }
};//updateDisplay;

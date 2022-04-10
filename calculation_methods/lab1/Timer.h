//
// Created by emely on 10.04.2022.
//

#ifndef LAB1__TIMER_H_
#define LAB1__TIMER_H_

#include <chrono>
class Timer {
 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> begin_;
 public:
  void start();
  long long int get_time_in_microseconds();
};

#endif //LAB1__TIMER_H_

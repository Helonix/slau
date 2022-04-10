//
// Created by emely on 10.04.2022.
//

#include "Timer.h"
void Timer::start() {
    begin_ = std::chrono::high_resolution_clock::now();
}
long long int Timer::get_time_in_microseconds() {
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(end - begin_).count();
}

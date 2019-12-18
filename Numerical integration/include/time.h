#pragma once

#include <iostream>
#include <chrono>

class Timer {
private:
  std::chrono::time_point<std::chrono::steady_clock> start, end;
  std::chrono::duration<double> duration;
public:
  Timer() {
    start = std::chrono::high_resolution_clock::now();
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
  }

  ~Timer() {
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Time: " << duration << std::endl;
  }
};
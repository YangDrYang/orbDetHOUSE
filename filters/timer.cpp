#include "timer.hpp"

Timer::Timer() {
    start = clock();
}

void Timer::tick() {
    start = clock();
}

double Timer::tock() {
    clock_t stop = clock();
    return ((double) (stop - start)) / CLOCKS_PER_SEC;
}

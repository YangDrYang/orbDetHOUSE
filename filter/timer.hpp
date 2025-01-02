#ifndef TIMER_H
#define TIMER_H

#include <ctime>

class Timer {

    private:

        clock_t start;

    public:

        Timer();

        void tick();

        double tock();

};

#endif

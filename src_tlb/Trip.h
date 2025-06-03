#ifndef TRIP
#define TRIP

#include <cstdint>

class Trip
{
public:
    Trip() : start_no(-1), end_no(-1), start_t(-1), end_t(-1), idx(-1) {}
    int16_t start_no;
    int16_t end_no;
    int32_t start_t;
    int32_t end_t;
    int idx;
    
    void Show()
    {
        printf("Trip start_no:%d end_no:%d start_t:%d end_t:%d idx:%d\n", start_no, end_no, start_t, end_t, idx);
    }
};


#endif
#include <unistd.h>
#include <stdio.h>
volatile char _lock_byte;

inline char TAS(char v) { return __sync_lock_test_and_set(&_lock_byte, v); } 
inline void SET(char v) { _lock_byte = v; }



inline void backoff(int *current_sleep_time_in_us){
    //printf("sleep for %d\n", current_sleep_time_in_us);
    usleep(*current_sleep_time_in_us);
    #ifdef DEF_SLEEP_MULTIPLIER
    *current_sleep_time_in_us = *current_sleep_time_in_us * DEF_SLEEP_MULTIPLIER;
    #else
    *current_sleep_time_in_us= *current_sleep_time_in_us ;//* 2;
    #endif
}

void lock() {
    #ifdef DEF_SLEEP_TIME 
    int current_sleep_time_in_us = DEF_SLEEP_TIME;
    #else
    int current_sleep_time_in_us = 2;
    #endif
	while (1) {
		while(_lock_byte);
		if (!TAS(1)) 
                {
                    return;
                }
		else backoff(&current_sleep_time_in_us);
	}
}

void unlock() {
	SET(0);
}



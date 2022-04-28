#if !defined TIMER_H
#define TIMER_H

#include <ctime>
#include <cmath>
//Returns the current time 
//Pre: None
//Post: Returns the current time.
time_t getTime ();

//Returns the difference between two variables of type time_h
//Pre: Start is a time variable that occurs before end
//Post: returns the amount of elapsed time between start and end.
double totalTime (time_t start, time_t end);

#endif
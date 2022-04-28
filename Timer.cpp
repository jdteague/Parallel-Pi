#include <chrono>
#include "Timer.h"

time_t getTime ()
{
	return time (NULL);
}

double totalTime (time_t start, time_t end)
{
	return difftime (end, start);
}
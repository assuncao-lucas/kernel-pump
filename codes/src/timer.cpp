#include "timer.h"

Timestamp::Timestamp()
{
}

Timestamp::~Timestamp()
{
}

Timer::Timer()
{
}

Timer::~Timer()
{
}

#ifdef _WIN32

WinTimestamp::WinTimestamp()
{
}

WinTimestamp::WinTimestamp(LARGE_INTEGER time)
{
	this->time_ = time;
}

WinTimestamp::~WinTimestamp()
{
}

LARGE_INTEGER WinTimestamp::time()
{
	return this->time_;
}

void WinTimestamp::set_time(LARGE_INTEGER time)
{
	this->time_ = time;
}

WinTimer::WinTimer()
{
	QueryPerformanceFrequency(&(this->frequency_));
}

WinTimer::~WinTimer()
{
}

void WinTimer::Clock(Timestamp* ts1)
{
    if(ts1 == NULL) return;
	LARGE_INTEGER t1;
	QueryPerformanceCounter(&t1);
	((WinTimestamp*)ts1)->set_time(t1);
}

double WinTimer::ElapsedTime(Timestamp* t1, Timestamp* t2)
{
	WinTimestamp * wt1 = (WinTimestamp*)t1;
	WinTimestamp * wt2 = (WinTimestamp*)t2;

	return (((wt2->time()).QuadPart - (wt1->time()).QuadPart)*1.0)/ (this->frequency_).QuadPart;
}

double WinTimer::CurrentElapsedTime(Timestamp* t1)
{
	WinTimestamp * wt1 = (WinTimestamp*)t1;
	WinTimestamp wt2;

	this->Clock(&wt2);

	return (((wt2.time()).QuadPart - (wt1->time()).QuadPart)*1.0)/ (this->frequency_).QuadPart;
}

#endif

#ifdef __unix__

#define BILLION 1E9

LinTimestamp::LinTimestamp()
{
}

LinTimestamp::LinTimestamp(struct timespec time)
{
	this->time_ = time;
}

LinTimestamp::~LinTimestamp()
{
}

void LinTimestamp::set_time(struct timespec time)
{
	this->time_ = time;
}

struct timespec LinTimestamp::time()
{
	return this->time_;
}

LinTimer::LinTimer()
{
}

LinTimer::~LinTimer()
{
}

void LinTimer::Clock(Timestamp* ts1)
{
    if(ts1 == nullptr) return;
	struct timespec t1;
	clock_gettime(CLOCK_REALTIME, &t1);
	((LinTimestamp*)ts1)->set_time(t1);
}

double LinTimer::ElapsedTime(Timestamp* ts1, Timestamp* ts2)
{
	struct timespec lt1 = ((LinTimestamp*)ts1)->time();
	struct timespec lt2 = ((LinTimestamp*)ts2)->time();

	return (lt2.tv_sec - lt1.tv_sec) + (lt2.tv_nsec - lt1.tv_nsec) / BILLION;
}

double LinTimer::CurrentElapsedTime(Timestamp* ts1)
{
    LinTimestamp ts2;
    this->Clock(&ts2);

	struct timespec lt1 = ((LinTimestamp*)ts1)->time();
	struct timespec lt2 = ts2.time();

	return (lt2.tv_sec - lt1.tv_sec) + (lt2.tv_nsec - lt1.tv_nsec) / BILLION;
}

#endif

static Timer * timer = NULL;

Timestamp * NewTimestamp(void)
{
#ifdef _WIN32
    return new WinTimestamp();
#endif

#ifdef __unix__
    return new LinTimestamp();
#endif
}

Timer * GetTimer(void)
{
    if(timer == NULL)
    {
#ifdef _WIN32
     timer = new WinTimer();
#endif
#ifdef __unix__
    timer = new LinTimer();
#endif
    }

return timer;

}

void DeleteTimer(void)
{
    if(timer != NULL)
    {
        delete timer;
        timer = NULL;
    }
}

#pragma once

#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __unix__
#include <sys/time.h>
#endif

class Timestamp
{
public:
	Timestamp();
	virtual ~Timestamp();
};

Timestamp *NewTimestamp(void);

class Timer
{
public:
	Timer();
	virtual ~Timer();
	virtual void Clock(Timestamp *) = 0;
	virtual double ElapsedTime(Timestamp *t1, Timestamp *t2) = 0;
	virtual double CurrentElapsedTime(Timestamp *t1) = 0;
	static Timer *timer_;
};

Timer *GetTimer(void);
void DeleteTimer(void);

#ifdef _WIN32

class WinTimestamp : public Timestamp
{
public:
	WinTimestamp();
	WinTimestamp(LARGE_INTEGER time);
	virtual ~WinTimestamp();
	LARGE_INTEGER time();
	void set_time(LARGE_INTEGER);

private:
	LARGE_INTEGER time_;
};

class WinTimer : public Timer
{
public:
	WinTimer();
	~WinTimer();
	virtual void Clock(Timestamp *);
	virtual double ElapsedTime(Timestamp *t1, Timestamp *t2);
	virtual double CurrentElapsedTime(Timestamp *t1);

private:
	LARGE_INTEGER frequency_;
};

#endif

#ifdef __unix__

class LinTimestamp : public Timestamp
{
public:
	LinTimestamp();
	LinTimestamp(struct timespec time);
	virtual ~LinTimestamp();
	struct timespec time();
	void set_time(struct timespec);

private:
	struct timespec time_;
};

class LinTimer : public Timer
{
public:
	LinTimer();
	~LinTimer();
	virtual void Clock(Timestamp *);
	virtual double ElapsedTime(Timestamp *t1, Timestamp *t2);
	virtual double CurrentElapsedTime(Timestamp *t1);
};
#endif

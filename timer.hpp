#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <string>

// Generic timer class wrapping C++11's std::chrono. Starts automatically when
// constructed; use one of the TimeElapsed functions to get readings. If the
// timer is stopped, each gives the elapsed time when it was stopped in the
// appropriate units; otherwise each gives the elapsed time at its invocation.
// Using Start() on an existing Timer object causes it to restart and makes the
// previous timing information unavailable.

class Timer {
	public:
		Timer();
		void Start();
		void Stop();
		long TimeElapsedInNanoseconds() const;
		long TimeElapsedInMicroseconds() const;
		long TimeElapsedInMilliseconds() const;
		long TimeElapsedInSeconds() const;
		long TimeElapsedInMinutes() const;
		long TimeElapsedInHours() const;
		std::string TimeElapsedInWords() const;

	private:
		std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
		std::chrono::time_point<std::chrono::high_resolution_clock> endTime;
		bool stopped;
		template<typename T> long TimeElapsedInUnit() const;
};

inline Timer::Timer() {
	Start();
}

inline void Timer::Start() {
	startTime = std::chrono::high_resolution_clock::now();
	stopped = false;
}

inline void Timer::Stop() {
	endTime = std::chrono::high_resolution_clock::now();
	stopped = true;
}

template<typename T>
inline long Timer::TimeElapsedInUnit() const {
	if(!stopped){
		//std::cerr << "Warning: asked for time elapsed in a timer which has not "
			//<< "been stopped." << std::endl;
		return std::chrono::duration_cast<T>(
				std::chrono::high_resolution_clock::now() - startTime).count();
	}
	return std::chrono::duration_cast<T>(endTime - startTime).count();
}

inline long Timer::TimeElapsedInNanoseconds() const {
	return TimeElapsedInUnit<std::chrono::nanoseconds>();
}

inline long Timer::TimeElapsedInMicroseconds() const {
	return TimeElapsedInUnit<std::chrono::microseconds>();
}

inline long Timer::TimeElapsedInMilliseconds() const {
	return TimeElapsedInUnit<std::chrono::milliseconds>();
}

inline long Timer::TimeElapsedInSeconds() const {
	return TimeElapsedInUnit<std::chrono::seconds>();
}

inline long Timer::TimeElapsedInMinutes() const {
	return TimeElapsedInUnit<std::chrono::minutes>();
}

inline long Timer::TimeElapsedInHours() const {
	return TimeElapsedInUnit<std::chrono::hours>();
}

inline std::string Timer::TimeElapsedInWords() const {
	long milliseconds = TimeElapsedInMilliseconds();

	int hours = milliseconds / 3600000;
	milliseconds = milliseconds % 3600000;

	int minutes = milliseconds / 60000;
	milliseconds = milliseconds % 60000;

	int seconds = milliseconds / 1000;
	milliseconds = milliseconds % 1000;

	std::string ret = "";
	if(hours > 0) ret += std::to_string(hours) + "h ";
	if(minutes > 0) ret += std::to_string(minutes) + "m ";
	if(seconds > 0) ret += std::to_string(seconds) + "s ";
	if(milliseconds > 0) ret += std::to_string(milliseconds) + "ms ";
	if(ret.empty()){
		ret = "0ms";
	} else {
		ret.erase(ret.end()-1);
	}
	return ret;
}
#endif

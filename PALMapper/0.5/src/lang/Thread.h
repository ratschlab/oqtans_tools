#pragma once

#include <pthread.h>
#include <unistd.h>

namespace lang {

class Mutex {

	friend class Signal;

public:

	class Locker {
	public:
		Locker(Mutex &mutex) : _mutex(mutex) {
			_mutex.lock();
		}

		~Locker() {
			_mutex.unlock();
		}

	private:
		Mutex &_mutex;
	};

	Mutex() {
		_mutex = _mutexInitializer;
	}

	void lock() {
		::pthread_mutex_lock(&_mutex);
	}

	void unlock() {
		::pthread_mutex_unlock(&_mutex);
	}

private:
	::pthread_mutex_t _mutex;
	static ::pthread_mutex_t _mutexInitializer;
};

class Signal {
public:
	Signal() {
		_cond = _initializer;
	}

	void wait(Mutex &mutex) {
		::pthread_cond_wait(&_cond, &mutex._mutex);
	}

	void notify() {
		::pthread_cond_signal(&_cond);
	}

	void notifyAll() {
		::pthread_cond_broadcast(&_cond);
	}

private:
	::pthread_cond_t _cond;
	static ::pthread_cond_t _initializer;
};

/**
 *
 */
class Runnable {
public:
	virtual void run() = 0;
};

/**
 *
 */
template <class T> class MethodRunnable : public Runnable {
public:
	MethodRunnable(void (T::*startMethod)()) {
		_startMethod = startMethod;
	}

	void run() {
		_startMethod();
	}

private:
	void (T::*_startMethod)();
};

/**
 *
 */
class Thread : Runnable {

public:

	Thread() {
		_runnable = this;
	}

	Thread(Runnable &runnable) {
		_runnable = &runnable;
	}

   virtual ~Thread() {
   }

   virtual void run() {
	   if (_runnable != NULL)
		   _runnable->run();
   }

   int launch();

   void join() {
	   ::pthread_join(_handle, NULL);
   }

   static void sleep(int milliseconds) {
	   ::sleep(milliseconds);
   }

   static void terminate();


private:
   Runnable *_runnable;
   static void *launchInterface(void *self);
   ::pthread_t _handle;
};

}

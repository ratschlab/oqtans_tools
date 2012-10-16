#include <stdio.h>
#include <new>
#include <lang/Thread.h>

namespace lang {

::pthread_mutex_t Mutex::_mutexInitializer = PTHREAD_MUTEX_INITIALIZER;
::pthread_cond_t Signal::_initializer = PTHREAD_COND_INITIALIZER;

int Thread::launch() {
   return ::pthread_create(&_handle, NULL, &launchInterface, this);
}

void *Thread::launchInterface(void *self) {
   try {
      ((Thread*) self)->run();
   } catch (std::exception &ex) {
	   ::printf("Unhandled exception in thread: %s", ex.what());
   } catch (...) {
      ::printf("Unhandled exception in thread\n");
   }
   return NULL;
}

}

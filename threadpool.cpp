#include "threadpool.h"

#include <unistd.h>
#include <iostream>

#define DISABLE_THREADS 0

Threadpool::RunningContext::RunningContext(Threadpool &threadpool)
    : _threadpool(threadpool)
{
    _threadpool._running = true;

    for (int i = 0; i < _threadpool._max_num_threads; i++) {
	_threadpool._threads.push_back(
	    new std::thread(
		std::bind(&Threadpool::_worker, &_threadpool)));
    }
}

Threadpool::RunningContext::~RunningContext()
{
  std::cout << "waiting for threads to finish\n";
    _threadpool._waitToFinish();
    std::cout << "done waiting\n";
}

Threadpool::Threadpool()
    : _max_num_threads(DISABLE_THREADS ?
		         0 :
   		         std::thread::hardware_concurrency() - 0),
      _num_running_threads_and_pending_tasks(0),
      _running(false)
{
}

void
Threadpool::scheduleOrExecute(const Task &task)
{
    if (_running) {
	std::scoped_lock lock(_mutex);
	if (_num_running_threads_and_pending_tasks < _max_num_threads) {
	    _num_running_threads_and_pending_tasks++;
	    _pending_tasks.push_back(task);
	    return;
	}
    }
    task();
}

void
Threadpool::_worker()
{
    while(true) {
	Task task;
	{
	    std::scoped_lock lock(_mutex);

	    if (_pending_tasks.empty()) {
		if (not _running) {
		    return;
		}
	    } else {
		task = _pending_tasks.back();
		_pending_tasks.pop_back();
	    }
	}

	if (task) {
	    task();
	    {
		std::scoped_lock lock(_mutex);
		_num_running_threads_and_pending_tasks--;
	    }
	} else {
	    usleep(10);
	}
    }
}

void
Threadpool::_waitToFinish()
{
    while (_running) {
	{
	    std::scoped_lock lock(_mutex);
	    if (_num_running_threads_and_pending_tasks == 0) {
		_running = false;
		break;
	    }
	}
	usleep(10000);
    }

    for (int i = 0; i < _max_num_threads; i++) {
	_threads[i]->join();
	delete _threads[i];
    }
    
    _threads.resize(0);
}

#ifndef __THREADPOOL_H__
#define __THREADPOOL_H__

#include <thread>
#include <mutex>
#include <functional>
#include <vector>

class Threadpool
{
public:
    Threadpool();

    typedef std::function<void ()> Task;
    void scheduleOrExecute(const Task &task);

    class RunningContext {
    public:
	RunningContext(Threadpool &threadpool);
	~RunningContext();

    private:
	Threadpool &_threadpool;
    };

private:
    void _worker();
    void _waitToFinish();

    std::vector<std::thread*> _threads;
    const int _max_num_threads;
    int _num_running_threads_and_pending_tasks;
    std::vector<Task> _pending_tasks;
    bool _running;
    
    std::mutex _mutex;
};

#endif

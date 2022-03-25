#ifndef __THREAD_SAFE_SET_H__
#define __THREAD_SAFE_SET_H__

#include <set>
#include <mutex>

template<class T> class ThreadSafeSet
{
public:
    const std::set<T> &get_set_unsafe() { return _set; };

    // For compatibility with set's insert.
    std::pair<int, bool> insert(const T &element) {
	std::scoped_lock lock(_mutex);
	return std::make_pair(0, _set.insert(element).second);
    }

private:
    std::set<T> _set;
    std::mutex _mutex;
};

#endif

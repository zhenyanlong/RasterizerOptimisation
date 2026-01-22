#pragma once

#include <atomic>
#include <thread>

class SpinLock {
private:
	std::atomic_flag flag = ATOMIC_FLAG_INIT;

public:
	void lock() {
		
		int backoff = 1;
		while (flag.test_and_set(std::memory_order_acquire)) {
			
			for (int i = 0; i < backoff && flag.test(std::memory_order_relaxed); ++i) {
				
#ifdef __x86_64__
				__builtin_ia32_pause();
#elif defined(_MSC_VER) && defined(_M_IX86)
				_mm_pause();
#endif
			}

			
			backoff = std::min(backoff * 2, 1024);

			
			if (backoff >= 64) {
				std::this_thread::yield();
			}
		}
	}

	void unlock() {
		flag.clear(std::memory_order_release);
	}

	
	bool try_lock() {
		return !flag.test_and_set(std::memory_order_acquire);
	}
};


class SpinLockGuard {
private:
	SpinLock& lock;

public:
	explicit SpinLockGuard(SpinLock& lock) : lock(lock) {
		lock.lock();
	}

	~SpinLockGuard() {
		lock.unlock();
	}

	SpinLockGuard(const SpinLockGuard&) = delete;
	SpinLockGuard& operator=(const SpinLockGuard&) = delete;
};
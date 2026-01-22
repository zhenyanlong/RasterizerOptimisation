#pragma once

#include <atomic>
#include <thread>

// 简单的自旋锁实现
class SpinLock {
private:
	std::atomic_flag flag = ATOMIC_FLAG_INIT;

public:
	void lock() {
		// 使用exponential backoff策略
		int backoff = 1;
		while (flag.test_and_set(std::memory_order_acquire)) {
			// 如果锁被占用，先自旋一段时间
			for (int i = 0; i < backoff && flag.test(std::memory_order_relaxed); ++i) {
				// 在x86上，使用pause指令减少功耗和总线争用
#ifdef __x86_64__
				__builtin_ia32_pause();
#elif defined(_MSC_VER) && defined(_M_IX86)
				_mm_pause();
#endif
			}

			// 指数退避，避免过多线程同时争抢
			backoff = std::min(backoff * 2, 1024);

			// 如果自旋一定次数后仍然无法获取锁，让出时间片
			if (backoff >= 64) {
				std::this_thread::yield();
			}
		}
	}

	void unlock() {
		flag.clear(std::memory_order_release);
	}

	// 尝试获取锁，非阻塞
	bool try_lock() {
		return !flag.test_and_set(std::memory_order_acquire);
	}
};

// 自旋锁的RAII包装器（类似std::lock_guard）
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

	// 禁止拷贝和移动
	SpinLockGuard(const SpinLockGuard&) = delete;
	SpinLockGuard& operator=(const SpinLockGuard&) = delete;
};
#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>
#include <stdexcept>

// 轻量级线程池类（适配渲染场景）
class ThreadPool {
private:
	// 工作线程列表
	std::vector<std::thread> workers;
	// 任务队列
	std::queue<std::function<void()>> tasks;

	// 同步相关
	std::mutex queue_mutex;       // 保护任务队列的锁
	std::condition_variable cv;   // 唤醒线程的条件变量
	bool stop;                    // 线程池停止标志
	std::atomic<int> task_count;  // 未完成的任务数（用于每帧等待）

public:
	// 构造函数：创建指定数量的工作线程
	explicit ThreadPool(size_t num_threads) : stop(false), task_count(0) {
		if (num_threads == 0) {
			num_threads = std::thread::hardware_concurrency(); // 默认使用CPU核心数
			if (num_threads == 0) num_threads = 4; // 兜底
		}

		// 创建工作线程
		for (size_t i = 0; i < num_threads; ++i) {
			workers.emplace_back([this]() {
				// 线程循环：等待并执行任务
				while (true) {
					std::function<void()> task;

					// 加锁获取任务
					{
						std::unique_lock<std::mutex> lock(this->queue_mutex);
						// 等待条件：有任务 或 线程池停止
						this->cv.wait(lock, [this]() {
							return this->stop || !this->tasks.empty();
							});

						// 线程池停止且无任务 → 退出
						if (this->stop && this->tasks.empty()) {
							return;
						}

						// 取出任务
						task = std::move(this->tasks.front());
						this->tasks.pop();
					}

					// 执行任务
					task();
					// 任务完成：计数减1
					this->task_count.fetch_sub(1);
				}
				});
		}
	}

	// 禁用拷贝构造和赋值
	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

	// 提交任务到线程池
	template<class F>
	void enqueue(F&& f) {
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			// 线程池已停止 → 拒绝提交任务
			if (stop) {
				throw std::runtime_error("enqueue on stopped ThreadPool");
			}
			// 添加任务到队列
			tasks.emplace(std::forward<F>(f));
			// 任务计数加1
			task_count.fetch_add(1);
		}
		// 唤醒一个等待的线程
		cv.notify_one();
	}

	// 等待所有已提交的任务执行完成（每帧渲染关键）
	void wait_all_tasks() {
		while (task_count.load() > 0) {
			std::this_thread::yield(); // 让出CPU，避免忙等
		}
	}

	// 销毁线程池：停止所有线程并等待退出
	~ThreadPool() {
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			stop = true;
		}
		// 唤醒所有线程
		cv.notify_all();
		// 等待所有线程退出
		for (std::thread& worker : workers) {
			if (worker.joinable()) {
				worker.join();
			}
		}
	}

	// 获取线程池的线程数量
	size_t get_thread_count() const {
		return workers.size();
	}
};
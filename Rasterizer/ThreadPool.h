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


class ThreadPool {
private:
	// List of worker threads
	std::vector<std::thread> workers;
	// Task queue
	std::queue<std::function<void()>> tasks;

	
	std::mutex queue_mutex;       // Lock for protecting the task queue
	std::condition_variable cv;   // The condition variable for waking up a thread
	bool stop;                    // Thread pool stop flag
	std::atomic<int> task_count;  // Number of unfinished tasks (used for per-frame waiting)

public:
	// Constructor: create a specified number of worker threads
	explicit ThreadPool(size_t num_threads) : stop(false), task_count(0) {
		if (num_threads == 0) {
			num_threads = std::thread::hardware_concurrency(); // Default to number of CPU cores
			if (num_threads == 0) num_threads = 4; // Fallback
		}

		// Create worker threads
		for (size_t i = 0; i < num_threads; ++i) {
			workers.emplace_back([this]() {
				// Thread loop: wait for and execute tasks
				while (true) {
					std::function<void()> task;

					// Lock and get task
					{
						std::unique_lock<std::mutex> lock(this->queue_mutex);
						// Wait condition: there are tasks or the thread pool is stopped
						this->cv.wait(lock, [this]() {
							return this->stop || !this->tasks.empty();
							});

						// Thread pool stopped and no tasks ¡ú exit
						if (this->stop && this->tasks.empty()) {
							return;
						}

						// Get task
						task = std::move(this->tasks.front());
						this->tasks.pop();
					}

					// Execute task
					task();
					// Task completed: decrement count
					this->task_count.fetch_sub(1);
				}
				});
		}
	}

	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

	// Submit a task to the thread pool
	template<class F>
	void enqueue(F&& f) {
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			// Thread pool stopped ¡ú reject task submission
			if (stop) {
				throw std::runtime_error("enqueue on stopped ThreadPool");
			}
			// Add task to queue	
			tasks.emplace(std::forward<F>(f));
			// Increment task count
			task_count.fetch_add(1);
		}
		// Wake up one waiting thread
		cv.notify_one();
	}

	// Wait for all submitted tasks to complete
	void wait_all_tasks() {
		while (task_count.load() > 0) {
			std::this_thread::yield(); // Yield CPU to avoid busy waiting
		}
	}

	
	~ThreadPool() {
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			stop = true;
		}
		
		cv.notify_all();
		
		for (std::thread& worker : workers) {
			if (worker.joinable()) {
				worker.join();
			}
		}
	}

	size_t get_thread_count() const {
		return workers.size();
	}
};
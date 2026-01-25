#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <atomic>
#include <condition_variable>
#include <memory>
#include <functional>
#include <future>
#include <deque>

// Simple support class for a 2D vector
class vec2D {
public:
	float x, y;

	// Default constructor initializes both components to 0
	vec2D() { x = y = 0.f; };

	// Constructor initializes components with given values
	vec2D(float _x, float _y) : x(_x), y(_y) {}

	// Constructor initializes components from a vec4
	vec2D(vec4 v) {
		x = v[0];
		y = v[1];
	}

	// Display the vector components
	void display() { std::cout << x << '\t' << y << std::endl; }

	// Overloaded subtraction operator for vector subtraction
	vec2D operator- (vec2D& v) {
		vec2D q;
		q.x = x - v.x;
		q.y = y - v.y;
		return q;
	}
};

// Class representing a triangle for rendering purposes
class triangle {
	Vertex v[3];       // Vertices of the triangle
	float area;        // Area of the triangle
	colour col[3];     // Colors for each vertex of the triangle

public:
	// Constructor initializes the triangle with three vertices
	triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
		v[0] = v1;
		v[1] = v2;
		v[2] = v3;

		// Calculate the 2D area of the triangle
		vec2D e1 = vec2D(v[1].p - v[0].p);
		vec2D e2 = vec2D(v[2].p - v[0].p);
		area = std::fabs(e1.x * e2.y - e1.y * e2.x);
	}

	// Helper function to compute the cross product for barycentric coordinates
	float getC(vec2D v1, vec2D v2, vec2D p) {
		vec2D e = v2 - v1;
		vec2D q = p - v1;
		return q.y * e.x - q.x * e.y;
	}

	// Compute barycentric coordinates for a given point
	bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
		alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
		beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
		gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;

		if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
		return true;
	}

	// Template function to interpolate values using barycentric coordinates
	template <typename T>
	T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
		return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
	}

	// Draw the triangle on the canvas
	void draw(Renderer& renderer, Light& L, float ka, float kd) {
		vec2D minV, maxV;
		getBoundsWindow(renderer.canvas, minV, maxV);

		if (area < 1.f) return;

		L.omega_i.normalise();

		// AVX optimization for 8 pixels at a time
		__m256 epsilon = _mm256_set1_ps(0.001f);
		__m256 zero = _mm256_setzero_ps();

		// Pre-compute coefficients for barycentric coordinate calculation
		float A_alpha = (v[0].p[1] - v[1].p[1]) / area;
		float B_alpha = (v[1].p[0] - v[0].p[0]) / area;
		float C_alpha = (v[0].p[0] * v[1].p[1] - v[1].p[0] * v[0].p[1]) / area;

		float A_beta = (v[1].p[1] - v[2].p[1]) / area;
		float B_beta = (v[2].p[0] - v[1].p[0]) / area;
		float C_beta = (v[1].p[0] * v[2].p[1] - v[2].p[0] * v[1].p[1]) / area;

		float A_gamma = (v[2].p[1] - v[0].p[1]) / area;
		float B_gamma = (v[0].p[0] - v[2].p[0]) / area;
		float C_gamma = (v[2].p[0] * v[0].p[1] - v[0].p[0] * v[2].p[1]) / area;

		int min_y = static_cast<int>(minV.y);
		int max_y = static_cast<int>(std::ceil(maxV.y));
		int min_x = static_cast<int>(minV.x);
		int max_x = static_cast<int>(std::ceil(maxV.x));
		int canvas_width = renderer.canvas.getWidth();

		// 使用互斥锁保护渲染操作
		static std::mutex renderMutex;

		// Process 8 pixels in parallel using SIMD
		for (int y = min_y; y <= max_y; y++) {
			float fy = static_cast<float>(y);

			for (int x = min_x; x <= max_x; x += 8) {
				int end_x = std::min(x + 7, canvas_width - 1);

				float xs[8];
				for (int i = 0; i < 8; i++) {
					int px = x + i;
					xs[i] = (px <= end_x) ? static_cast<float>(px) : static_cast<float>(x);
				}

				__m256 x_vec = _mm256_loadu_ps(xs);
				__m256 y_vec = _mm256_set1_ps(fy);

				// Compute barycentric coordinates using SIMD
				__m256 alpha = _mm256_add_ps(
					_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(A_alpha), x_vec),
						_mm256_mul_ps(_mm256_set1_ps(B_alpha), y_vec)),
					_mm256_set1_ps(C_alpha)
				);
				__m256 beta = _mm256_add_ps(
					_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(A_beta), x_vec),
						_mm256_mul_ps(_mm256_set1_ps(B_beta), y_vec)),
					_mm256_set1_ps(C_beta)
				);
				__m256 gamma = _mm256_add_ps(
					_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(A_gamma), x_vec),
						_mm256_mul_ps(_mm256_set1_ps(B_gamma), y_vec)),
					_mm256_set1_ps(C_gamma)
				);

				// Create a mask for pixels inside the triangle
				__m256 alpha_ge_zero = _mm256_cmp_ps(alpha, zero, _CMP_GE_OS);
				__m256 beta_ge_zero = _mm256_cmp_ps(beta, zero, _CMP_GE_OS);
				__m256 gamma_ge_zero = _mm256_cmp_ps(gamma, zero, _CMP_GE_OS);
				__m256 mask_inside = _mm256_and_ps(
					_mm256_and_ps(alpha_ge_zero, beta_ge_zero),
					gamma_ge_zero
				);

				if (_mm256_movemask_ps(mask_inside) == 0) continue;

				// Interpolate depth using barycentric coordinates
				__m256 depth = _mm256_add_ps(
					_mm256_add_ps(_mm256_mul_ps(alpha, _mm256_set1_ps(v[0].p[2])),
						_mm256_mul_ps(beta, _mm256_set1_ps(v[1].p[2]))),
					_mm256_mul_ps(gamma, _mm256_set1_ps(v[2].p[2]))
				);

				// Load Z-buffer values for the 8 pixels
				float zbuf_vals[8];
				for (int i = 0; i < 8; i++) {
					int px = x + i;
					zbuf_vals[i] = (px <= end_x) ? renderer.zbuffer(px, y) : 1e9f;
				}
				__m256 zbuf = _mm256_loadu_ps(zbuf_vals);

				// Perform depth test using SIMD
				__m256 depth_test1 = _mm256_cmp_ps(zbuf, depth, _CMP_GT_OS);
				__m256 depth_test2 = _mm256_cmp_ps(depth, epsilon, _CMP_GT_OS);
				__m256 depth_test = _mm256_and_ps(depth_test1, depth_test2);
				__m256 final_mask = _mm256_and_ps(mask_inside, depth_test);
				int mask = _mm256_movemask_ps(final_mask);

				// Store barycentric coordinates and depth to arrays
				float alpha_arr[8], beta_arr[8], gamma_arr[8];
				_mm256_storeu_ps(alpha_arr, alpha);
				_mm256_storeu_ps(beta_arr, beta);
				_mm256_storeu_ps(gamma_arr, gamma);
				float depth_arr[8];
				_mm256_storeu_ps(depth_arr, depth);

				// 保护像素写入操作
				std::lock_guard<std::mutex> lock(renderMutex);

				// Process each pixel based on the final mask
				for (int i = 0; i < 8; i++) {
					if (mask & (1 << i)) {
						int px = x + i;
						if (px > end_x) break;

						// Interpolate color and normal for the pixel
						colour c = interpolate(alpha_arr[i], beta_arr[i], gamma_arr[i],
							v[0].rgb, v[1].rgb, v[2].rgb);
						c.clampColour();
						vec4 normal = interpolate(alpha_arr[i], beta_arr[i], gamma_arr[i],
							v[0].normal, v[1].normal, v[2].normal);
						normal.normalise();

						// Compute lighting using the interpolated normal
						float dot_val = std::max(vec4::dot(L.omega_i, normal), 0.0f);
						colour diffuse = (c * kd) * (L.L * dot_val);
						colour ambient = L.ambient * ka;
						colour final_col = diffuse + ambient;

						// Convert final color to RGB and draw the pixel
						unsigned char r, g, b;
						final_col.toRGB(r, g, b);
						renderer.canvas.draw(px, y, r, g, b);
						renderer.zbuffer(px, y) = depth_arr[i];
					}
				}
			}
		}
	}

	// Draw only within a specific tile (用于多线程)
	void drawTile(Renderer& renderer, Light& L, float ka, float kd,
		int tileX, int tileY, int tileWidth, int tileHeight) {
		vec2D minV, maxV;
		getBoundsWindow(renderer.canvas, minV, maxV);

		// Clip to tile bounds
		minV.x = std::max(minV.x, static_cast<float>(tileX));
		minV.y = std::max(minV.y, static_cast<float>(tileY));
		maxV.x = std::min(maxV.x, static_cast<float>(tileX + tileWidth));
		maxV.y = std::min(maxV.y, static_cast<float>(tileY + tileHeight));

		if (minV.x >= maxV.x || minV.y >= maxV.y) return;
		if (area < 1.f) return;

		// 直接调用draw函数，它已经有边界裁剪
		draw(renderer, L, ka, kd);
	}

	// Compute the 2D bounds of the triangle
	void getBounds(vec2D& minV, vec2D& maxV) {
		minV = vec2D(v[0].p);
		maxV = vec2D(v[0].p);
		for (unsigned int i = 1; i < 3; i++) {
			minV.x = std::min(minV.x, v[i].p[0]);
			minV.y = std::min(minV.y, v[i].p[1]);
			maxV.x = std::max(maxV.x, v[i].p[0]);
			maxV.y = std::max(maxV.y, v[i].p[1]);
		}
	}

	// Compute the 2D bounds of the triangle, clipped to the canvas
	void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
		getBounds(minV, maxV);
		minV.x = std::max(minV.x, static_cast<float>(0));
		minV.y = std::max(minV.y, static_cast<float>(0));
		maxV.x = std::min(maxV.x, static_cast<float>(canvas.getWidth()));
		maxV.y = std::min(maxV.y, static_cast<float>(canvas.getHeight()));
	}

	void drawBounds(GamesEngineeringBase::Window& canvas) {
		vec2D minV, maxV;
		getBounds(minV, maxV);

		for (int y = (int)minV.y; y < (int)maxV.y; y++) {
			for (int x = (int)minV.x; x < (int)maxV.x; x++) {
				canvas.draw(x, y, 255, 0, 0);
			}
		}
	}

	void display() {
		for (unsigned int i = 0; i < 3; i++) {
			v[i].p.display();
		}
		std::cout << std::endl;
	}
};

// ==================== 简化但高效的线程池实现 ====================

// 简化的工作队列
class WorkQueue {
private:
	typedef std::function<void()> TaskType;
	std::deque<TaskType> tasks;
	mutable std::mutex mutex;

public:
	WorkQueue() = default;

	void push(TaskType task) {
		std::lock_guard<std::mutex> lock(mutex);
		tasks.push_back(std::move(task));
	}

	bool try_pop(TaskType& task) {
		std::lock_guard<std::mutex> lock(mutex);
		if (tasks.empty()) {
			return false;
		}
		task = std::move(tasks.front());
		tasks.pop_front();
		return true;
	}

	bool empty() const {
		std::lock_guard<std::mutex> lock(mutex);
		return tasks.empty();
	}
};

// 高性能线程池（简化版，避免C++17特性）
class ThreadPool {
private:
	std::atomic_bool done;
	WorkQueue queue;
	std::vector<std::thread> threads;
	std::condition_variable condition;
	mutable std::mutex queue_mutex;

public:
	ThreadPool() : done(false) {
		unsigned const thread_count = std::thread::hardware_concurrency();
		try {
			for (unsigned i = 0; i < thread_count; ++i) {
				threads.push_back(std::thread(&ThreadPool::worker, this));
			}
		}
		catch (...) {
			done = true;
			throw;
		}
	}

	~ThreadPool() {
		done = true;
		condition.notify_all();
		for (auto& thread : threads) {
			if (thread.joinable()) {
				thread.join();
			}
		}
	}

	template<typename FunctionType>
	void submit(FunctionType f) {
		{
			std::lock_guard<std::mutex> lock(queue_mutex);
			queue.push(std::function<void()>(f));
		}
		condition.notify_one();
	}

	void wait_all() {
		while (!queue.empty()) {
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}

		// 等待所有任务完成
		bool has_tasks = true;
		while (has_tasks) {
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
			std::lock_guard<std::mutex> lock(queue_mutex);
			has_tasks = !queue.empty();
		}
	}

private:
	void worker() {
		while (!done) {
			std::function<void()> task;
			{
				std::unique_lock<std::mutex> lock(queue_mutex);
				condition.wait(lock, [this] {
					return done || !queue.empty();
					});

				if (done && queue.empty()) {
					return;
				}

				if (!queue.try_pop(task)) {
					continue;
				}
			}

			if (task) {
				task();
			}
		}
	}
};

// 高性能多线程渲染器
class HighPerformanceRenderer {
private:
	// 静态线程池实例（单例模式）
	static ThreadPool& getThreadPool() {
		static ThreadPool pool;
		return pool;
	}

public:
	// 高性能多线程渲染
	static void renderMultiThread(Renderer& renderer, const std::vector<triangle>& triangles,
		Light& L, float ka, float kd, int tileSize = 32) {

		if (triangles.empty()) return;

		ThreadPool& pool = getThreadPool();
		int canvasWidth = renderer.canvas.getWidth();
		int canvasHeight = renderer.canvas.getHeight();

		// 计算tile网格
		int tilesX = (canvasWidth + tileSize - 1) / tileSize;
		int tilesY = (canvasHeight + tileSize - 1) / tileSize;

		// 为每个三角形和每个tile创建任务
		for (const auto& tri : triangles) {
			// 获取三角形的边界
			vec2D minV, maxV;
			const_cast<triangle&>(tri).getBoundsWindow(renderer.canvas, minV, maxV);

			if (minV.x >= maxV.x || minV.y >= maxV.y) continue;

			// 计算三角形覆盖的tile范围
			int startTileX = std::max(0, static_cast<int>(minV.x) / tileSize);
			int endTileX = std::min(tilesX - 1, static_cast<int>(maxV.x) / tileSize);
			int startTileY = std::max(0, static_cast<int>(minV.y) / tileSize);
			int endTileY = std::min(tilesY - 1, static_cast<int>(maxV.y) / tileSize);

			// 为三角形覆盖的每个tile创建渲染任务
			for (int ty = startTileY; ty <= endTileY; ++ty) {
				for (int tx = startTileX; tx <= endTileX; ++tx) {
					int tileX = tx * tileSize;
					int tileY = ty * tileSize;
					int tileWidth = std::min(tileSize, canvasWidth - tx * tileSize);
					int tileHeight = std::min(tileSize, canvasHeight - ty * tileSize);

					// 提交任务到线程池（使用lambda捕获所需参数）
					pool.submit([tri, &renderer, &L, ka, kd, tileX, tileY, tileWidth, tileHeight]() {
						// 创建三角形副本以避免const问题
						triangle triCopy = tri;
						triCopy.drawTile(renderer, L, ka, kd, tileX, tileY, tileWidth, tileHeight);
						});
				}
			}
		}

		// 等待所有任务完成
		pool.wait_all();
	}

	// 按物体划分的渲染（减少任务数量）
	static void renderByObject(Renderer& renderer, const std::vector<triangle>& triangles,
		Light& L, float ka, float kd) {

		if (triangles.empty()) return;

		ThreadPool& pool = getThreadPool();

		// 每个任务处理一个三角形
		for (const auto& tri : triangles) {
			pool.submit([tri, &renderer, &L, ka, kd]() {
				triangle triCopy = tri;
				triCopy.draw(renderer, L, ka, kd);
				});
		}

		// 等待所有任务完成
		pool.wait_all();
	}

	// 批量渲染（减少任务提交开销）
	static void renderBatch(Renderer& renderer, const std::vector<triangle>& triangles,
		Light& L, float ka, float kd, int batchSize = 100) {

		if (triangles.empty()) return;

		ThreadPool& pool = getThreadPool();

		// 分批处理三角形
		int totalTriangles = static_cast<int>(triangles.size());

		for (int i = 0; i < totalTriangles; i += batchSize) {
			int endIdx = std::min(i + batchSize, totalTriangles);

			pool.submit([&renderer, &triangles, &L, ka, kd, i, endIdx]() {
				for (int j = i; j < endIdx; ++j) {
					triangle triCopy = triangles[j];
					triCopy.draw(renderer, L, ka, kd);
				}
				});
		}

		// 等待所有任务完成
		pool.wait_all();
	}
};
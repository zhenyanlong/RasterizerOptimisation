#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>
#include <atomic>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"
#include "ThreadPool.h"

// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.
void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Combine perspective, camera, and world transformations for the mesh
    matrix p = renderer.perspective * camera * mesh->world;


    // Iterate through all triangles in the mesh
    for (triIndices& ind : mesh->triangles) {
        Vertex t[3]; // Temporary array to store transformed triangle vertices

        // Transform each vertex of the triangle
        #pragma loop(3)
        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
            t[i].p.divideW(); // Perspective division to normalize coordinates

            // Transform normals into world space for accurate lighting
            // no need for perspective correction as no shearing or non-uniform scaling
            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal; 
            t[i].normal.normalise();

            // Map normalized device coordinates to screen space
            t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
            t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
            t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

            // Copy vertex colours
            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        // Clip triangles with Z-values outside [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        // Create a triangle object and render it
        triangle tri(t[0], t[1], t[2]);
        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

std::atomic<int> mesh_index(0);
std::atomic<int> tri_index(0);
ThreadPool thread_pool(22); // 传0自动适配CPU核心数

// another type of render function, transmit a vector array of meshes to render
void renderScene(Renderer& renderer, std::vector<Mesh*>& scene, matrix& camera, Light& L, int num_threads) {
	// Combine perspective, camera, and world transformations for the mesh
	// create threads
	// 先收集所有需要处理的三角形任务
	struct RenderTask {
		Mesh* mesh;
		size_t triangle_index;
		matrix transform;
	};

	std::vector<RenderTask> tasks;

	// 预计算所有任务
	for (auto& mesh : scene) {
		if (!mesh) continue;

		matrix p = renderer.perspective * camera * mesh->world;

		for (size_t i = 0; i < mesh->triangles.size(); ++i) {
			tasks.push_back({ mesh, i, p });
		}
	}

	// 使用原子索引来分配任务
	std::atomic<size_t> task_index(0);
	size_t total_tasks = tasks.size();

	auto render_task = [&](int thread_id) {
		while (true) {
			// 原子获取下一个任务索引
			size_t current_task = task_index.fetch_add(1, std::memory_order_relaxed);

			if (current_task >= total_tasks) {
				break;
			}

			const auto& task = tasks[current_task];
			triIndices& ind = task.mesh->triangles[task.triangle_index];
			Vertex t[3];

			// Transform vertices
#pragma loop(3)
			for (unsigned int i = 0; i < 3; i++) {
				t[i].p = task.transform * task.mesh->vertices[ind.v[i]].p;
				t[i].p.divideW();

				t[i].normal = task.mesh->world * task.mesh->vertices[ind.v[i]].normal;
				t[i].normal.normalise();

				t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
				t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
				t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1];

				t[i].rgb = task.mesh->vertices[ind.v[i]].rgb;
			}

			// Clip
			if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

			// Draw
			triangle tri(t[0], t[1], t[2]);
			tri.draw(renderer, L, task.mesh->ka, task.mesh->kd);
		}
		};

	// 创建线程
	std::vector<std::thread> threads;
	threads.reserve(num_threads);

	for (int i = 0; i < num_threads; ++i) {
		threads.emplace_back(render_task, i);
	}

	// 等待完成
	for (auto& thread : threads) {
		if (thread.joinable()) {
			thread.join();
		}
	}
	

	// 线程执行函数：每个线程循环获取网格，处理该网格的所有三角形
	//auto render_worker = [&]() {
	//	while (true) {
	//		// 1. 原子获取当前要处理的网格索引（递增前的值）
	//		int mesh_idx = current_mesh_index.fetch_add(1);
	//		// 2. 检查是否超出场景范围，超出则退出线程
	//		if (mesh_idx >= static_cast<int>(scene.size())) {
	//			break;
	//		}
	//		// 3. 获取当前网格（确保不会越界）
	//		Mesh* mesh = scene[mesh_idx];
	//		if (!mesh) continue; // 防御性检查：空指针跳过

	//		// 4. 计算该网格的变换矩阵（每个网格只计算一次，避免重复计算）
	//		matrix p = renderer.perspective * camera * mesh->world;

	//		// 5. 处理该网格的所有三角形
	//		for (triIndices& ind : mesh->triangles) {
	//			Vertex t[3];
 //               #pragma unroll(3)
	//			for (unsigned int i = 0; i < 3; i++) {
	//				t[i].p = p * mesh->vertices[ind.v[i]].p;
	//				t[i].p.divideW();
	//				t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
	//				t[i].normal.normalise();
	//				t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
	//				t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
	//				t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1];
	//				t[i].rgb = mesh->vertices[ind.v[i]].rgb;
	//			}
	//			// 裁剪Z值超出范围的三角形
	//			if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

	//			// 渲染三角形（Renderer的锁会保护资源）
	//			triangle tri(t[0], t[1], t[2]);
	//			tri.draw(renderer, L, mesh->ka, mesh->kd);
	//		}
	//	}
	//	};

	//// 创建并启动线程
	//std::vector<std::thread> threads;
	//threads.reserve(multi_nums);
	//for (int i = 0; i < multi_nums; i++) {
	//	threads.emplace_back(render_worker);
	//}

	//// 等待所有线程执行完成（关键：必须join，否则渲染未完成就继续）
	//for (auto& t : threads) {
	//	if (t.joinable()) {
	//		t.join();
	//	}
	//}
}
void renderScene(Renderer& renderer, std::vector<Mesh*>& scene, matrix& camera, Light& L, ThreadPool& thread_pool) {
	// 遍历所有网格，为每个网格创建渲染任务并提交到线程池
	for (size_t mesh_idx = 0; mesh_idx < scene.size(); ++mesh_idx) {
		Mesh* mesh = scene[mesh_idx];
		if (!mesh) continue; // 空指针防御

		// 提交渲染单个网格的任务（捕获值，避免引用失效）
		thread_pool.enqueue([&renderer, mesh, &camera, &L]() {
			// 计算该网格的变换矩阵（每个网格仅计算一次）
			matrix p = renderer.perspective * camera * mesh->world;

			// 处理该网格的所有三角形
			for (triIndices& ind : mesh->triangles) {
				Vertex t[3];
#pragma loop(3)
				for (unsigned int i = 0; i < 3; i++) {
					t[i].p = p * mesh->vertices[ind.v[i]].p;
					t[i].p.divideW();
					t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
					t[i].normal.normalise();
					// 屏幕空间映射
					t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
					t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
					t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1];
					t[i].rgb = mesh->vertices[ind.v[i]].rgb;
				}

				// 裁剪Z值超出范围的三角形
				if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) {
					continue;
				}

				// 渲染三角形（Renderer的全局锁保护资源）
				triangle tri(t[0], t[1], t[2]);
				tri.draw(renderer, L, mesh->ka, mesh->kd);
			}
			});
	}

	// 等待当前帧所有渲染任务完成（关键：确保渲染完再present）
	thread_pool.wait_all_tasks();
}

// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
    Renderer renderer;
    // create light source {direction, diffuse intensity, ambient intensity}
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };
    // camera is just a matrix
    matrix camera = matrix::makeIdentity(); // Initialize the camera with identity matrix

    bool running = true; // Main loop control variable

    std::vector<Mesh*> scene; // Vector to store scene objects

    // Create a sphere and a rectangle mesh
    Mesh mesh = Mesh::makeSphere(1.0f, 10, 20);
    //Mesh mesh2 = Mesh::makeRectangle(-2, -1, 2, 1);

    // add meshes to scene
    scene.push_back(&mesh);
   // scene.push_back(&mesh2); 

    float x = 0.0f, y = 0.0f, z = -4.0f; // Initial translation parameters
    mesh.world = matrix::makeTranslation(x, y, z);
    //mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput(); // Handle user input
        renderer.clear(); // Clear the canvas for the next frame

        // Apply transformations to the meshes
     //   mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
        mesh.world = matrix::makeTranslation(x, y, z);

        // Handle user inputs for transformations
        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
        if (renderer.canvas.keyPressed('A')) x += -0.1f;
        if (renderer.canvas.keyPressed('D')) x += 0.1f;
        if (renderer.canvas.keyPressed('W')) y += 0.1f;
        if (renderer.canvas.keyPressed('S')) y += -0.1f;
        if (renderer.canvas.keyPressed('Q')) z += 0.1f;
        if (renderer.canvas.keyPressed('E')) z += -0.1f;

        // Render each object in the scene
        for (auto& m : scene)
            render(renderer, m, camera, L);

        renderer.present(); // Display the rendered frame
    }
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
    unsigned int r = rng.getRandomInt(0, 3);

    switch (r) {
    case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
    default: return matrix::makeIdentity();
    }
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    bool running = true;

    std::vector<Mesh*> scene;

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        /*for (auto& m : scene)
            render(renderer, m, camera, L);*/
		renderScene(renderer, scene, camera, L, thread_pool);
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    std::vector<Mesh*> scene;

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.push_back(m);
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.push_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        // Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		/*for (auto& m : scene)
			render(renderer, m, camera, L);*/
		renderScene(renderer, scene, camera, L, thread_pool);
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene: 10x10x10 cubes + random rotation each cube + camera circular orbiting
// No input variables
void scene3() {
	Renderer renderer;
	matrix camera = matrix::makeIdentity();
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

	std::vector<Mesh*> scene;
	struct rRot { float x; float y; float z; };
	std::vector<rRot> cubeRotations;

	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

	// Create a 10x10x10 cubic grid (totaling 1000 cells)
	const int gridSize = 10;          // Grid dimensions
	const float cubeSize = 0.3f;      // Cube size (to avoid excessive overlap)
	const float cubePadding = 0.4f; // Cube spacing
	const float gridCenter = (gridSize - 1) * (cubeSize + cubePadding) * 0.5f; // Grid center offset
	const float cameraDistance = 8.0f; // Camera distance to grid center
	for (int z = 0; z < gridSize; z++) {
		for (int y = 0; y < gridSize; y++) {
			for (int x = 0; x < gridSize; x++) {

				Mesh* cube = new Mesh();
				*cube = Mesh::makeCube(cubeSize);
				scene.push_back(cube);

				float posX = (x * (cubeSize + cubePadding)) - gridCenter;
				float posY = (y * (cubeSize + cubePadding)) - gridCenter;
				float posZ = (z * (cubeSize + cubePadding)) - gridCenter;
				cube->world = matrix::makeTranslation(posX, posY, posZ);

				// Generate random rotation increments for each cube
				rRot rot{
					rng.getRandomFloat(-0.02f, 0.02f),  // X-axis rotation increment
					rng.getRandomFloat(-0.02f, 0.02f),  // Y-axis rotation increment
					rng.getRandomFloat(-0.02f, 0.02f)   // Z-axis rotation increment
				};
				cubeRotations.push_back(rot);

				// Initial random rotation
				cube->world = cube->world * matrix::makeRotateXYZ(
					rng.getRandomFloat(0.f, 3.14f),
					rng.getRandomFloat(0.f, 3.14f),
					rng.getRandomFloat(0.f, 3.14f)
				);
			}
		}
	}

	// Camera rotation control variables
	float cameraAngle = 0.0f;         // Camera rotation angle (radians)
	const float angleStep = 0.01f;    // Rotation increment per frame (controls rotation speed)
	const float fullCircle = 2 * 3.1415926f; // One full circle in radians (360 degrees)

	
	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0; // Cycle count

	bool running = true;
	while (running) {
		renderer.canvas.checkInput();
		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		renderer.clear();


		// Update cube rotations
		for (unsigned int i = 0; i < cubeRotations.size(); i++) {
			scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(
				cubeRotations[i].x,
				cubeRotations[i].y,
				cubeRotations[i].z
			);
		}

		// Update camera position (circular orbit around Y-axis)    
		// Calculate camera's circular coordinates
		float camX = sin(cameraAngle) * cameraDistance;
        float camZ = cos(cameraAngle) * cameraDistance; 
        
		// Construct camera matrix
		matrix camRotate = matrix::makeRotateY(-cameraAngle); 
		matrix camTranslate = matrix::makeTranslation(-camX, 0.0f, -(camZ));
		camera = camRotate*camTranslate;

		// Update camera rotation angle and check if a full circle is completed
		cameraAngle += angleStep;
		if (cameraAngle >= fullCircle) {
			cameraAngle -= fullCircle; 
			cycle++;                   

			
			end = std::chrono::high_resolution_clock::now();
			double cycleTimeMs = std::chrono::duration<double, std::milli>(end - start).count();
			std::cout << cycle << " :" << cycleTimeMs << "ms\n";
			start = std::chrono::high_resolution_clock::now(); 
		}

		/*for (auto& m : scene) {
			render(renderer, m, camera, L);
		}*/

		renderScene(renderer, scene, camera, L, thread_pool);

		renderer.present();
	}
	for (auto& m : scene) {
		delete m;
	}
}

// Entry point of the application
// No input variables
int main() {
	
	//std::cout << "线程池创建完成，线程数：" << thread_pool.get_thread_count() << std::endl;
    // Uncomment the desired scene function to run
    //scene1();
    scene2();
    //scene3();
    //sceneTest(); 
    

    return 0;
}
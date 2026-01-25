#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"

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
        #pragma unroll(3)
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

// 多线程渲染函数 - 使用高性能线程池
void renderMultiThread(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L,
	int mode = 0) {  // 0: tile-based, 1: by object, 2: batch
	matrix p = renderer.perspective * camera * mesh->world;

	std::vector<triangle> triangles;

	for (triIndices& ind : mesh->triangles) {
		Vertex t[3];

		for (unsigned int i = 0; i < 3; i++) {
			t[i].p = p * mesh->vertices[ind.v[i]].p;
			t[i].p.divideW();

			t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
			t[i].normal.normalise();

			t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
			t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
			t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1];

			t[i].rgb = mesh->vertices[ind.v[i]].rgb;
		}

		if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

		triangles.emplace_back(t[0], t[1], t[2]);
	}

	// 根据模式选择不同的渲染策略
	switch (mode) {
	case 0:  // Tile-based（最适合大量小三角形）
		HighPerformanceRenderer::renderMultiThread(renderer, triangles, L, mesh->ka, mesh->kd, 32);
		break;
	case 1:  // By object（减少任务数量）
		HighPerformanceRenderer::renderByObject(renderer, triangles, L, mesh->ka, mesh->kd);
		break;
	case 2:  // Batch（平衡任务开销）
		HighPerformanceRenderer::renderBatch(renderer, triangles, L, mesh->ka, mesh->kd, 50);
		break;
	default:
		// 默认使用batch模式
		HighPerformanceRenderer::renderBatch(renderer, triangles, L, mesh->ka, mesh->kd, 50);
		break;
	}
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

        for (auto& m : scene)
            render(renderer, m, camera, L);
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

        for (auto& m : scene)
            render(renderer, m, camera, L);
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
    int renderMode = 1;  // 默认使用tile-based多线程
	
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
		for (auto& m : scene) {
			if (renderMode == 0) {
				render(renderer, m, camera, L);
			}
			else {
				renderMultiThread(renderer, m, camera, L, renderMode - 1);
			}
		}

		renderer.present();
	}
	for (auto& m : scene) {
		delete m;
	}
}

// Entry point of the application
// No input variables
int main() {
    // Uncomment the desired scene function to run
    //scene1();
    //scene2();
    scene3();
    //sceneTest(); 
    

    return 0;
}
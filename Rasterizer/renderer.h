#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "GamesEngineeringBase.h"
#include "zbuffer.h"
#include "matrix.h"

struct Tile {
	int index;                    // Tile index
	int x, y;                     // Tile position (top-left corner)
	int width, height;            // Tile size
	std::vector<int> triangleIndices; // Indices of triangles belonging to this Tile

    Tile(int _index, int _x, int _y, int _w, int _h) : index(_index), x(_x), y(_y), width(_w), height(_h)
    {

    }

};

class TileManager {
private:
	int screenWidth, screenHeight;
	int tileSize;
	int numTilesX, numTilesY;
	std::vector<Tile> tiles;
	std::vector<triangle&> tileTriangles; 

public:
	void init(int width, int height, int tileSize = 32) {
		this->screenWidth = width;
		this->screenHeight = height;
		this->tileSize = tileSize;

		numTilesX = (width + tileSize - 1) / tileSize;
		numTilesY = (height + tileSize - 1) / tileSize;

		tiles.clear();
		tileTriangles.clear();
		for (int y = 0; y < numTilesY; ++y) {
			for (int x = 0; x < numTilesX; ++x) {
				int tileW = std::min(tileSize, width - x * tileSize-1);
				int tileH = std::min(tileSize, height - y * tileSize-1);
				
				tiles.emplace_back(y * numTilesX + x, x * tileSize, y * tileSize, tileW, tileH);
				
			}
		}
	}

	void clearTileTriangles() {
		tileTriangles.clear();
	}

	void addTriangleToTiles(const vec2D& minV, const vec2D& maxV, triangle& tri) {
		
		tileTriangles.emplace_back(tri);
		int triIndex = static_cast<int>(tileTriangles.size()) - 1; 
		
		int minTileX = std::max(0, static_cast<int>(minV.x) / tileSize);
		int minTileY = std::max(0, static_cast<int>(minV.y) / tileSize);
		int maxTileX = std::min(numTilesX - 1, static_cast<int>(maxV.x) / tileSize);
		int maxTileY = std::min(numTilesY - 1, static_cast<int>(maxV.y) / tileSize);

		
		for (int ty = minTileY; ty <= maxTileY; ++ty) {
			for (int tx = minTileX; tx <= maxTileX; ++tx) {
				int tileIndex = ty * numTilesX + tx;
				
				tiles[tileIndex].triangleIndices.push_back(triIndex);

			}
		}
	}
	
	int getTileCount() const { return tiles.size(); }
	Tile& getTile(int index) { return tiles[index]; }
	const std::vector<int>& getTileTriangles(int tileIndex) const {
		return tiles[tileIndex].triangleIndices;
	}

	void clear() {
		for (auto& tile : tiles) {
			tile.triangleIndices.clear();
		}
		tileTriangles.clear();
	}
};

// The `Renderer` class handles rendering operations, including managing the
// Z-buffer, canvas, and perspective transformations for a 3D scene.
class Renderer {
    float fov = 90.0f * M_PI / 180.0f; // Field of view in radians (converted from degrees)
    float aspect = 4.0f / 3.0f;        // Aspect ratio of the canvas (width/height)
    float n = 0.1f;                    // Near clipping plane distance
    float f = 100.0f;                  // Far clipping plane distance
public:
    Zbuffer<float> zbuffer;                  // Z-buffer for depth management
    GamesEngineeringBase::Window canvas;     // Canvas for rendering the scene
    matrix perspective;                      // Perspective projection matrix  
	TileManager tileManager;

    // Constructor initializes the canvas, Z-buffer, and perspective projection matrix.
    Renderer() {
        canvas.create(1024, 768, "Raster");  // Create a canvas with specified dimensions and title
        zbuffer.create(1024, 768);           // Initialize the Z-buffer with the same dimensions
        perspective = matrix::makePerspective(fov, aspect, n, f); // Set up the perspective matrix

		tileManager.init(1024, 768,256); // Initialize the tile manager with the canvas dimensions
    }

    // Clears the canvas and resets the Z-buffer.
    void clear() {
        canvas.clear();  // Clear the canvas (sets all pixels to the background color)
        zbuffer.clear(); // Reset the Z-buffer to the farthest depth
		tileManager.clear(); // Clear the tile manager's triangle lists
    }

    // Presents the current canvas frame to the display.
    void present() {
        canvas.present(); // Display the rendered frame
    }
};

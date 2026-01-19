#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <algorithm>
#include <cmath>

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
    // Input Variables:
    // - v1, v2, v3: Vertices defining the triangle
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
    // Input Variables:
    // - v1, v2: Edges defining the vector
    // - p: Point for which coordinates are being calculated
    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }

    // Compute barycentric coordinates for a given point
    // Input Variables:
    // - p: Point to check within the triangle
    // Output Variables:
    // - alpha, beta, gamma: Barycentric coordinates of the point
    // Returns true if the point is inside the triangle, false otherwise
    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;

        if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
        return true;
    }

    // Template function to interpolate values using barycentric coordinates
    // Input Variables:
    // - alpha, beta, gamma: Barycentric coordinates
    // - a1, a2, a3: Values to interpolate
    // Returns the interpolated value
    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }

    // Draw the triangle on the canvas
    // Input Variables:
    // - renderer: Renderer object for drawing
    // - L: Light object for shading calculations
    // - ka, kd: Ambient and diffuse lighting coefficients
    void draw(Renderer& renderer, Light& L, float ka, float kd) {
  //      vec2D minV, maxV;

  //      // Get the screen-space bounds of the triangle
  //      getBoundsWindow(renderer.canvas, minV, maxV);

  //      // Skip very small triangles
  //      if (area < 1.f) return;

		//// pre - normalise light direction
  //      L.omega_i.normalise();

  //      // Iterate over the bounding box and check each pixel
  //      for (int y = (int)(minV.y); y < (int)ceil(maxV.y); y++) {
  //          
  //          for (int x = (int)(minV.x); x < (int)ceil(maxV.x); x++) {
  //              float alpha, beta, gamma;

  //              // Check if the pixel lies inside the triangle
  //              if (getCoordinates(vec2D((float)x, (float)y), alpha, beta, gamma)) {
  //                  // Interpolate color, depth, and normals
  //                  colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
  //                  c.clampColour();
  //                  float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
  //                  vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
  //                  normal.normalise();

  //                  // Perform Z-buffer test and apply shading
  //                  if (renderer.zbuffer(x, y) > depth && depth > 0.001f) {
  //                      // typical shader begin
  //                      //L.omega_i.normalise();
  //                      float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
  //                      colour a = (c * kd) * (L.L * dot) + (L.ambient * ka); // using kd instead of ka for ambient
  //                      // typical shader end
  //                      unsigned char r, g, b;
  //                      a.toRGB(r, g, b);
  //                      renderer.canvas.draw(x, y, r, g, b);
  //                      renderer.zbuffer(x, y) = depth;
  //                  }
  //              }
  //               
  //              
  //          }
  //      }
        // -- SIMD 4 pixel -- //
  //      vec2D minV, maxV;
  //      getBoundsWindow(renderer.canvas, minV, maxV);

  //  
  //      if (area < 1.f) return;

  //      L.omega_i.normalise();

  //      __m128 kd_vec = _mm_set1_ps(kd);
  //      __m128 ka_vec = _mm_set1_ps(ka);
  //      __m128 epsilon = _mm_set1_ps(0.001f);
  //      __m128 zero = _mm_setzero_ps();

		//// Precompute coefficients for barycentric coordinate calculation
  //      float A_alpha = (v[0].p[1] - v[1].p[1]) / area;
  //      float B_alpha = (v[1].p[0] - v[0].p[0]) / area;
  //      float C_alpha = (v[0].p[0] * v[1].p[1] - v[1].p[0] * v[0].p[1]) / area;

	 //   float A_beta = (v[1].p[1] - v[2].p[1]) / area;
	 //   float B_beta = (v[2].p[0] - v[1].p[0]) / area;
  //      float C_beta = (v[1].p[0] * v[2].p[1] - v[2].p[0] * v[1].p[1]) / area;

  //      float A_gamma = (v[2].p[1] - v[0].p[1]) / area;
  //      float B_gamma = (v[0].p[0] - v[2].p[0]) / area;
  //      float C_gamma = (v[2].p[0] * v[0].p[1] - v[0].p[0] * v[2].p[1]) / area;

		//// Loop over the bounding box of the triangle
  //      int min_y = static_cast<int>(minV.y);
  //      int max_y = static_cast<int>(std::ceil(maxV.y));
  //      int min_x = static_cast<int>(minV.x);
  //      int max_x = static_cast<int>(std::ceil(maxV.x));
  //      int canvas_width = renderer.canvas.getWidth();

		//// Process 4 pixels in parallel using SIMD
  //      for (int y = min_y; y <= max_y; y++) {
  //          float fy = static_cast<float>(y);

  //          for (int x = min_x; x <= max_x; x += 4) {
  //              
  //              int end_x = std::min(x + 3, canvas_width - 1);
  //              
  //              float xs[4];
  //              xs[0] = static_cast<float>(x);
  //              xs[1] = (x+1 <= end_x) ? static_cast<float>(x+1) : static_cast<float>(x);
  //              xs[2] = (x+2 <= end_x) ? static_cast<float>(x+2) : static_cast<float>(x);
  //              xs[3] = (x+3 <= end_x) ? static_cast<float>(x+3) : static_cast<float>(x);

  //              __m128 x_vec = _mm_load_ps(xs);
  //              __m128 y_vec = _mm_set1_ps(fy);

		//		// Compute barycentric coordinates using SIMD
  //              __m128 alpha = _mm_add_ps(
  //                  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(A_alpha), x_vec), _mm_mul_ps(_mm_set1_ps(B_alpha), y_vec)),
  //                  _mm_set1_ps(C_alpha)
  //              );
  //              __m128 beta = _mm_add_ps(
  //                  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(A_beta), x_vec), _mm_mul_ps(_mm_set1_ps(B_beta), y_vec)),
  //                  _mm_set1_ps(C_beta)
  //              );
  //              __m128 gamma = _mm_add_ps(
  //                  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(A_gamma), x_vec), _mm_mul_ps(_mm_set1_ps(B_gamma), y_vec)),
  //                  _mm_set1_ps(C_gamma)
  //              );

		//		// Create a mask for pixels inside the triangle
  //              __m128 mask_inside = _mm_and_ps(
  //                  _mm_and_ps(_mm_cmpge_ps(alpha, zero), _mm_cmpge_ps(beta, zero)),
  //                  _mm_cmpge_ps(gamma, zero)
  //              );

  //              if (_mm_movemask_ps(mask_inside) == 0) continue;

		//		// Interpolate depth using barycentric coordinates
  //              __m128 depth = _mm_add_ps(
  //                  _mm_add_ps(_mm_mul_ps(alpha, _mm_set1_ps(v[0].p[2])),
  //                             _mm_mul_ps(beta, _mm_set1_ps(v[1].p[2]))),
  //                  _mm_mul_ps(gamma, _mm_set1_ps(v[2].p[2]))
  //              );

		//		// Load Z-buffer values for the 4 pixels
  //              float zbuf_vals[4];
  //              for (int i=0; i<4; i++) {
  //                  int px = x + i;
  //                  zbuf_vals[i] = (px <= end_x) ? renderer.zbuffer(px, y) : 1e9f;
  //              }
  //              __m128 zbuf = _mm_load_ps(zbuf_vals);

		//		// Perform depth test using SIMD
  //              __m128 depth_test1 = _mm_cmpgt_ps(zbuf, depth);
  //              __m128 depth_test2 = _mm_cmpgt_ps(depth, epsilon);
  //              __m128 depth_test = _mm_and_ps(depth_test1, depth_test2);
  //              __m128 final_mask = _mm_and_ps(mask_inside, depth_test);
  //              int mask = _mm_movemask_ps(final_mask);

		//		// Store barycentric coordinates and depth to arrays for per-pixel processing
  //              float alpha_arr[4], beta_arr[4], gamma_arr[4];
  //              _mm_store_ps(alpha_arr, alpha);
  //              _mm_store_ps(beta_arr, beta);
  //              _mm_store_ps(gamma_arr, gamma);
  //              float depth_arr[4];
  //              _mm_store_ps(depth_arr, depth);

		//		// Process each pixel based on the final mask
  //              #pragma unroll(4)
  //              for (int i = 0; i < 4; i++) {
  //                  if (mask & (1 << i)) {
  //                      int px = x + i;
  //                      if (px > end_x) break;

		//				// Interpolate color and normal for the pixel
  //                      colour c = interpolate(alpha_arr[i], beta_arr[i], gamma_arr[i], v[0].rgb, v[1].rgb, v[2].rgb);
  //                      c.clampColour();
  //                      vec4 normal = interpolate(alpha_arr[i], beta_arr[i], gamma_arr[i], v[0].normal, v[1].normal, v[2].normal);
  //                      normal.normalise();

		//				// Compute lighting using the interpolated normal
  //                      float dot_val = std::max(vec4::dot(L.omega_i, normal), 0.0f);
  //                      colour diffuse = (c * kd) * (L.L * dot_val);
  //                      colour ambient = L.ambient * ka;
  //                      colour final_col = diffuse + ambient;

		//				// Convert final color to RGB and draw the pixel
  //                      unsigned char r, g, b;
  //                      final_col.toRGB(r, g, b);
  //                      renderer.canvas.draw(px, y, r, g, b);
  //                      renderer.zbuffer(px, y) = depth_arr[i];
  //                  }
  //              }
  //          }
  //      }
        vec2D minV, maxV;
        getBoundsWindow(renderer.canvas, minV, maxV);
        
        if (area < 1.f) return;
        
        L.omega_i.normalise();
        
        
        __m256 kd_vec = _mm256_set1_ps(kd);
        __m256 ka_vec = _mm256_set1_ps(ka);
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
        
        // Loop over the bounding box of the triangle
        int min_y = static_cast<int>(minV.y);
        int max_y = static_cast<int>(std::ceil(maxV.y));
        int min_x = static_cast<int>(minV.x);
        int max_x = static_cast<int>(std::ceil(maxV.x));
        int canvas_width = renderer.canvas.getWidth();
        
        // Process 8 pixels in parallel using SIMD
        for (int y = min_y; y <= max_y; y++) {
        	float fy = static_cast<float>(y);
        
        	// x步长从4改为8
        	for (int x = min_x; x <= max_x; x += 8) {
        
        		// 结束x从x+3改为x+7
        		int end_x = std::min(x + 7, canvas_width - 1);
        
        		// 数组大小从4改为8
        		float xs[8];
                #pragma unroll(8)
        		for (int i = 0; i < 8; i++) {
        			int px = x + i;
        			xs[i] = (px <= end_x) ? static_cast<float>(px) : static_cast<float>(x);
        		}
        
        		// 替换为AVX指令
        		__m256 x_vec = _mm256_loadu_ps(xs);
        		__m256 y_vec = _mm256_set1_ps(fy);
        
        		// Compute barycentric coordinates using SIMD
        		__m256 alpha = _mm256_add_ps(
        			_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(A_alpha), x_vec), _mm256_mul_ps(_mm256_set1_ps(B_alpha), y_vec)),
        			_mm256_set1_ps(C_alpha)
        		);
        		__m256 beta = _mm256_add_ps(
        			_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(A_beta), x_vec), _mm256_mul_ps(_mm256_set1_ps(B_beta), y_vec)),
        			_mm256_set1_ps(C_beta)
        		);
        		__m256 gamma = _mm256_add_ps(
        			_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(A_gamma), x_vec), _mm256_mul_ps(_mm256_set1_ps(B_gamma), y_vec)),
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
        
        		// 替换为_mm256_movemask_ps
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
        
        		// Process each pixel based on the final mask
                #pragma unroll(8) 
        		for (int i = 0; i < 8; i++) {
        			if (mask & (1 << i)) {
        				int px = x + i;
        				if (px > end_x) break;
        
        				// Interpolate color and normal for the pixel
        				colour c = interpolate(alpha_arr[i], beta_arr[i], gamma_arr[i], v[0].rgb, v[1].rgb, v[2].rgb);
        				c.clampColour();
        				vec4 normal = interpolate(alpha_arr[i], beta_arr[i], gamma_arr[i], v[0].normal, v[1].normal, v[2].normal);
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

    // Compute the 2D bounds of the triangle
    // Output Variables:
    // - minV, maxV: Minimum and maximum bounds in 2D space
    void getBounds(vec2D& minV, vec2D& maxV) {
        minV = vec2D(v[0].p);
        maxV = vec2D(v[0].p);
        #pragma unroll(2)
        for (unsigned int i = 1; i < 3; i++) {
            minV.x = std::min(minV.x, v[i].p[0]);
            minV.y = std::min(minV.y, v[i].p[1]);
            maxV.x = std::max(maxV.x, v[i].p[0]);
            maxV.y = std::max(maxV.y, v[i].p[1]);
        }
    }

    // Compute the 2D bounds of the triangle, clipped to the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    // Output Variables:
    // - minV, maxV: Clipped minimum and maximum bounds
    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        getBounds(minV, maxV);
        minV.x = std::max(minV.x, static_cast<float>(0));
        minV.y = std::max(minV.y, static_cast<float>(0));
        maxV.x = std::min(maxV.x, static_cast<float>(canvas.getWidth()));
        maxV.y = std::min(maxV.y, static_cast<float>(canvas.getHeight()));
    }

    // Debugging utility to display the triangle bounds on the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
                canvas.draw(x, y, 255, 0, 0);
            }
        }
    }

    // Debugging utility to display the coordinates of the triangle vertices
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }
};

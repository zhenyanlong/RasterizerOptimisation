#pragma once

#include <iostream>
#include <vector>
#include "vec4.h"

// Matrix class for 4x4 transformation matrices
class matrix {
    union {
        float m[4][4]; // 2D array representation of the matrix
        float a[16];   // 1D array representation of the matrix for linear access
    };

public:
    // Default constructor initializes the matrix as an identity matrix
    matrix() {
        identity();
    }

    // Access matrix elements by row and column
    float& operator()(unsigned int row, unsigned int col) { return m[row][col]; }

    // Display the matrix elements in a readable format
    void display() {
        for (unsigned int i = 0; i < 4; i++) {
            for (unsigned int j = 0; j < 4; j++)
                std::cout << m[i][j] << '\t';
            std::cout << std::endl;
        }
    }

    // Multiply the matrix by a 4D vector
    // Input Variables:
    // - v: vec4 object to multiply with the matrix
    // Returns the resulting transformed vec4
    vec4 operator * (const vec4& v) const {
        vec4 result;
        result[0] = a[0] * v[0] + a[1] * v[1] + a[2] * v[2] + a[3] * v[3];
        result[1] = a[4] * v[0] + a[5] * v[1] + a[6] * v[2] + a[7] * v[3];
        result[2] = a[8] * v[0] + a[9] * v[1] + a[10] * v[2] + a[11] * v[3];
        result[3] = a[12] * v[0] + a[13] * v[1] + a[14] * v[2] + a[15] * v[3];
        return result;
    }

    // Multiply the matrix by another matrix
    // Input Variables:
    // - mx: Another matrix to multiply with
    // Returns the resulting matrix
    matrix operator * (const matrix& mx) const {
        // roolling 
		/*matrix ret;
		for (int row = 0; row < 4; ++row) {
			for (int col = 0; col < 4; ++col) {
				ret.a[row * 4 + col] =
					a[row * 4 + 0] * mx.a[0 * 4 + col] +
					a[row * 4 + 1] * mx.a[1 * 4 + col] +
					a[row * 4 + 2] * mx.a[2 * 4 + col] +
					a[row * 4 + 3] * mx.a[3 * 4 + col];
			}
		}
		return ret;*/
        // UnRolling
		matrix ret;
		// 展开循环，使用局部变量减少内存访问
		// 第一行
		const float a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3];
		const float b00 = mx.a[0], b01 = mx.a[1], b02 = mx.a[2], b03 = mx.a[3];
		const float b10 = mx.a[4], b11 = mx.a[5], b12 = mx.a[6], b13 = mx.a[7];
		const float b20 = mx.a[8], b21 = mx.a[9], b22 = mx.a[10], b23 = mx.a[11];
		const float b30 = mx.a[12], b31 = mx.a[13], b32 = mx.a[14], b33 = mx.a[15];

		// 计算第一行
		ret.a[0] = a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
		ret.a[1] = a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
		ret.a[2] = a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
		ret.a[3] = a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

		// 第二行
		const float a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7];
		ret.a[4] = a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
		ret.a[5] = a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
		ret.a[6] = a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
		ret.a[7] = a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

		// 第三行
		const float a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11];
		ret.a[8] = a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
		ret.a[9] = a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
		ret.a[10] = a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
		ret.a[11] = a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

		// 第四行
		const float a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];
		ret.a[12] = a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
		ret.a[13] = a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
		ret.a[14] = a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
		ret.a[15] = a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;

		return ret;
    }

    // Create a perspective projection matrix
    // Input Variables:
    // - fov: Field of view in radians
    // - aspect: Aspect ratio of the viewport
    // - n: Near clipping plane
    // - f: Far clipping plane
    // Returns the perspective matrix
    static matrix makePerspective(float fov, float aspect, float n, float f) {
        // The last solution
		/*matrix m;
		m.zero();
		float tanHalfFov = std::tan(fov / 2.0f);

		m.a[0] = 1.0f / (aspect * tanHalfFov);
		m.a[5] = 1.0f / tanHalfFov;
		m.a[10] = -f / (f - n);
		m.a[11] = -(f * n) / (f - n);
		m.a[14] = -1.0f;
		return m;*/
        // Now
		matrix m;
		m.zero();

		// 预计算常用值，避免重复计算
		const float tanHalfFov = std::tan(fov * 0.5f);  // 提前计算半角tan值
		const float invTanHalfFov = 1.0f / tanHalfFov;  // 计算倒数，避免除法
		const float invAspectTanHalfFov = invTanHalfFov / aspect;
		const float range = n - f;
		const float invRange = 1.0f / range;  // 预计算倒数

		m.a[0] = invAspectTanHalfFov;
		m.a[5] = invTanHalfFov;
		m.a[10] = f * invRange;          // 优化：使用预计算的倒数
		m.a[11] = f * n * invRange;      // 优化：减少除法操作
		m.a[14] = -1.0f;

		return m;
    }

    // Create a translation matrix
    // Input Variables:
    // - tx, ty, tz: Translation amounts along the X, Y, and Z axes
    // Returns the translation matrix
    static matrix makeTranslation(float tx, float ty, float tz) {
        matrix m;
        m.identity();
        m.a[3] = tx;
        m.a[7] = ty;
        m.a[11] = tz;
        return m;
    }

    // Create a rotation matrix around the Z-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateZ(float aRad) {
		/*matrix m;
		m.identity();
		m.a[0] = std::cos(aRad);
		m.a[1] = -std::sin(aRad);
		m.a[4] = std::sin(aRad);
		m.a[5] = std::cos(aRad);
		return m;*/

        // Now
		matrix m;
		m.identity();
		const float c = std::cos(aRad);
		const float s = std::sin(aRad);
		m.a[0] = c;
		m.a[1] = -s;
		m.a[4] = s;
		m.a[5] = c;
		return m;
    }

    // Create a rotation matrix around the X-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateX(float aRad) {
		/*matrix m;
		m.identity();
		m.a[5] = std::cos(aRad);
		m.a[6] = -std::sin(aRad);
		m.a[9] = std::sin(aRad);
		m.a[10] = std::cos(aRad);
		return m;*/

		// Now
		matrix m;
		m.identity();
		const float c = std::cos(aRad);
		const float s = std::sin(aRad);
		m.a[5] = c;
		m.a[6] = -s;
		m.a[9] = s;
		m.a[10] = c;
		return m;
    }

    // Create a rotation matrix around the Y-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateY(float aRad) {
		/*matrix m;
		m.identity();
		m.a[0] = std::cos(aRad);
		m.a[2] = std::sin(aRad);
		m.a[8] = -std::sin(aRad);
		m.a[10] = std::cos(aRad);
		return m;*/

		matrix m;
		m.identity();
		const float c = std::cos(aRad);
		const float s = std::sin(aRad);
		m.a[0] = c;
		m.a[2] = s;
		m.a[8] = -s;
		m.a[10] = c;
		return m;
    }

    // Create a composite rotation matrix from X, Y, and Z rotations
    // Input Variables:
    // - x, y, z: Rotation angles in radians around each axis
    // Returns the composite rotation matrix
    static matrix makeRotateXYZ(float x, float y, float z) {
        //return matrix::makeRotateX(x) * matrix::makeRotateY(y) * matrix::makeRotateZ(z);

		// Now
		// 预计算三角函数值，避免重复计算
		const float cx = std::cos(x), sx = std::sin(x);
		const float cy = std::cos(y), sy = std::sin(y);
		const float cz = std::cos(z), sz = std::sin(z);

		matrix m;
		m.identity();

		// 组合旋转矩阵 (Rx * Ry * Rz)
		m.a[0] = cy * cz;
		m.a[1] = -cy * sz;
		m.a[2] = sy;

		m.a[4] = sx * sy * cz + cx * sz;
		m.a[5] = -sx * sy * sz + cx * cz;
		m.a[6] = -sx * cy;

		m.a[8] = -cx * sy * cz + sx * sz;
		m.a[9] = cx * sy * sz + sx * cz;
		m.a[10] = cx * cy;

		return m;
    }

    // Create a scaling matrix
    // Input Variables:
    // - s: Scaling factor
    // Returns the scaling matrix
    static matrix makeScale(float s) {
        matrix m;
        s = std::max(s, 0.01f); // Ensure scaling factor is not too small
        m.identity();
        m.a[0] = s;
        m.a[5] = s;
        m.a[10] = s;
        return m;
    }

    // Create an identity matrix
    // Returns an identity matrix
    static matrix makeIdentity() {
		/*matrix m;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				m.m[i][j] = (i == j) ? 1.0f : 0.0f;
			}
		}
		return m;*/

		// Now
		matrix m;
		// 展开循环
		m.a[0] = m.a[5] = m.a[10] = m.a[15] = 1.0f;
		m.a[1] = m.a[2] = m.a[3] = m.a[4] =
			m.a[6] = m.a[7] = m.a[8] = m.a[9] =
			m.a[11] = m.a[12] = m.a[13] = m.a[14] = 0.0f;
		return m;
    }

private:
    // Set all elements of the matrix to 0
    void zero() {
		/*for (unsigned int i = 0; i < 16; i++)
			a[i] = 0.f;*/

		// Now
		a[0] = a[1] = a[2] = a[3] =
			a[4] = a[5] = a[6] = a[7] =
			a[8] = a[9] = a[10] = a[11] =
			a[12] = a[13] = a[14] = a[15] = 0.0f;
    }

    // Set the matrix as an identity matrix
    void identity() {
		/*for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				m[i][j] = (i == j) ? 1.0f : 0.0f;
			}
		}*/

		// Now
		// 展开初始化，避免循环
		a[0] = a[5] = a[10] = a[15] = 1.0f;
		a[1] = a[2] = a[3] = a[4] =
			a[6] = a[7] = a[8] = a[9] =
			a[11] = a[12] = a[13] = a[14] = 0.0f;
    }
};



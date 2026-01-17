#pragma once

#include <iostream>
#include <vector>
#include "vec4.h"
#include <immintrin.h>  
#include <intrin.h>     

// Matrix class for 4x4 transformation matrices
class matrix {
    union {
		alignas(32)   float m[4][4]; // 2D array representation of the matrix
		alignas(32)   float a[16];   // 1D array representation of the matrix for linear access
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
    vec4 operator * (vec4& v) const {
		/*vec4 result;
		result[0] = a[0] * v[0] + a[1] * v[1] + a[2] * v[2] + a[3] * v[3];
		result[1] = a[4] * v[0] + a[5] * v[1] + a[6] * v[2] + a[7] * v[3];
		result[2] = a[8] * v[0] + a[9] * v[1] + a[10] * v[2] + a[11] * v[3];
		result[3] = a[12] * v[0] + a[13] * v[1] + a[14] * v[2] + a[15] * v[3];
		return result;*/


		//-- SIMD --//
		//// 将向量的每个分量广播到单独的__m128寄存器
		//__m128 vx = _mm_set1_ps(v[0]);
		//__m128 vy = _mm_set1_ps(v[1]);
		//__m128 vz = _mm_set1_ps(v[2]);
		//__m128 vw = _mm_set1_ps(v[3]);

		//__m128 row0 = _mm_set_ps(a[12], a[8], a[4], a[0]);
		//__m128 row1 = _mm_set_ps(a[13], a[9], a[5], a[1]);
		//__m128 row2 = _mm_set_ps(a[14], a[10], a[6], a[2]);
		//__m128 row3 = _mm_set_ps(a[15], a[11], a[7], a[3]);
		//__m128 result = _mm_add_ps(_mm_add_ps(_mm_mul_ps(row0, vx),
		//	_mm_mul_ps(row1, vy)),
		//	_mm_add_ps(_mm_mul_ps(row2, vz),
		//		_mm_mul_ps(row3, vw)));

		//// 存储结果
		//alignas(16) float result_arr[4];
		//_mm_store_ps(result_arr, result);

		//return vec4(result_arr[0], result_arr[1], result_arr[2], result_arr[3]);

		// -- SIMD 优化版 -- //
		// 将向量两两合并进行广播
		__m256 vxy = _mm256_set_m128(_mm_set1_ps(v[1]), _mm_set1_ps(v[0])); // v.y | v.x
		__m256 vzw = _mm256_set_m128(_mm_set1_ps(v[3]), _mm_set1_ps(v[2])); // v.w | v.z
		// 将矩阵的行两两合并
		__m256 row01 = _mm256_set_m128(_mm_set_ps(a[13], a[9], a[5], a[1]), _mm_set_ps(a[12], a[8], a[4], a[0])); // row1 | row0
		__m256 row23 = _mm256_set_m128(_mm_set_ps(a[15], a[11], a[7], a[3]), _mm_set_ps(a[14], a[10], a[6], a[2])); // row3 | row2
		// 将v和row进行乘法和加法运算
		__m256 mul0123 = _mm256_add_ps(_mm256_mul_ps(row01, vxy), _mm256_mul_ps(row23, vzw)); // 将结果相加
		// 将__256结果拆分回__m128以便存储
		__m128 result = _mm_add_ps(_mm256_castps256_ps128(mul0123), _mm256_extractf128_ps(mul0123, 1));

		// 储存结果
		alignas(16) float result_arr[4];
		_mm_store_ps(result_arr, result);
		return vec4(result_arr[0], result_arr[1], result_arr[2], result_arr[3]);

		
	
    }

    // Multiply the matrix by another matrix
    // Input Variables:
    // - mx: Another matrix to multiply with
    // Returns the resulting matrix
    matrix operator * (const matrix& mx) const {
        // rolling 
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
		//matrix ret;
		//// 展开循环，使用局部变量减少内存访问
		//// 第一行
		//const float a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3];
		//const float b00 = mx.a[0], b01 = mx.a[1], b02 = mx.a[2], b03 = mx.a[3];
		//const float b10 = mx.a[4], b11 = mx.a[5], b12 = mx.a[6], b13 = mx.a[7];
		//const float b20 = mx.a[8], b21 = mx.a[9], b22 = mx.a[10], b23 = mx.a[11];
		//const float b30 = mx.a[12], b31 = mx.a[13], b32 = mx.a[14], b33 = mx.a[15];

		//// 计算第一行
		//ret.a[0] = a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
		//ret.a[1] = a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
		//ret.a[2] = a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
		//ret.a[3] = a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

		//// 第二行
		//const float a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7];
		//ret.a[4] = a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
		//ret.a[5] = a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
		//ret.a[6] = a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
		//ret.a[7] = a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

		//// 第三行
		//const float a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11];
		//ret.a[8] = a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
		//ret.a[9] = a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
		//ret.a[10] = a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
		//ret.a[11] = a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

		//// 第四行
		//const float a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];
		//ret.a[12] = a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
		//ret.a[13] = a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
		//ret.a[14] = a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
		//ret.a[15] = a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;

		//return ret;



		


		// -- SIMD --//
		//matrix ret;

		//// 加载右矩阵mx的4行（矩阵乘法计算时需要的维度）
		//__m128 mx_row0 = _mm_load_ps(&mx.a[0]);   // mx[0][0], mx[0][1], mx[0][2], mx[0][3]
		//__m128 mx_row1 = _mm_load_ps(&mx.a[4]);   // mx[1][0], mx[1][1], mx[1][2], mx[1][3]
		//__m128 mx_row2 = _mm_load_ps(&mx.a[8]);   // mx[2][0], mx[2][1], mx[2][2], mx[2][3]
		//__m128 mx_row3 = _mm_load_ps(&mx.a[12]);  // mx[3][0], mx[3][1], mx[3][2], mx[3][3]

		//// ========== 循环展开：处理第0行（原col=0） ==========
		//__m128 a0_0 = _mm_set1_ps(a[0 * 4 + 0]);  // 广播左矩阵第0行第0列元素
		//__m128 a0_1 = _mm_set1_ps(a[0 * 4 + 1]);  // 广播左矩阵第0行第1列元素
		//__m128 a0_2 = _mm_set1_ps(a[0 * 4 + 2]);  // 广播左矩阵第0行第2列元素
		//__m128 a0_3 = _mm_set1_ps(a[0 * 4 + 3]);  // 广播左矩阵第0行第3列元素
		//__m128 result0 = _mm_add_ps(
		//	_mm_add_ps(_mm_mul_ps(a0_0, mx_row0), _mm_mul_ps(a0_1, mx_row1)),
		//	_mm_add_ps(_mm_mul_ps(a0_2, mx_row2), _mm_mul_ps(a0_3, mx_row3))
		//);
		//_mm_store_ps(&ret.a[0 * 4], result0);     // 存储结果矩阵第0行

		//// ========== 循环展开：处理第1行（原col=1） ==========
		//__m128 a1_0 = _mm_set1_ps(a[1 * 4 + 0]);  // 广播左矩阵第1行第0列元素
		//__m128 a1_1 = _mm_set1_ps(a[1 * 4 + 1]);  // 广播左矩阵第1行第1列元素
		//__m128 a1_2 = _mm_set1_ps(a[1 * 4 + 2]);  // 广播左矩阵第1行第2列元素
		//__m128 a1_3 = _mm_set1_ps(a[1 * 4 + 3]);  // 广播左矩阵第1行第3列元素
		//__m128 result1 = _mm_add_ps(
		//	_mm_add_ps(_mm_mul_ps(a1_0, mx_row0), _mm_mul_ps(a1_1, mx_row1)),
		//	_mm_add_ps(_mm_mul_ps(a1_2, mx_row2), _mm_mul_ps(a1_3, mx_row3))
		//);
		//_mm_store_ps(&ret.a[1 * 4], result1);     // 存储结果矩阵第1行

		//// ========== 循环展开：处理第2行（原col=2） ==========
		//__m128 a2_0 = _mm_set1_ps(a[2 * 4 + 0]);  // 广播左矩阵第2行第0列元素
		//__m128 a2_1 = _mm_set1_ps(a[2 * 4 + 1]);  // 广播左矩阵第2行第1列元素
		//__m128 a2_2 = _mm_set1_ps(a[2 * 4 + 2]);  // 广播左矩阵第2行第2列元素
		//__m128 a2_3 = _mm_set1_ps(a[2 * 4 + 3]);  // 广播左矩阵第2行第3列元素
		//__m128 result2 = _mm_add_ps(
		//	_mm_add_ps(_mm_mul_ps(a2_0, mx_row0), _mm_mul_ps(a2_1, mx_row1)),
		//	_mm_add_ps(_mm_mul_ps(a2_2, mx_row2), _mm_mul_ps(a2_3, mx_row3))
		//);
		//_mm_store_ps(&ret.a[2 * 4], result2);     // 存储结果矩阵第2行

		//// ========== 循环展开：处理第3行（原col=3） ==========
		//__m128 a3_0 = _mm_set1_ps(a[3 * 4 + 0]);  // 广播左矩阵第3行第0列元素
		//__m128 a3_1 = _mm_set1_ps(a[3 * 4 + 1]);  // 广播左矩阵第3行第1列元素
		//__m128 a3_2 = _mm_set1_ps(a[3 * 4 + 2]);  // 广播左矩阵第3行第2列元素
		//__m128 a3_3 = _mm_set1_ps(a[3 * 4 + 3]);  // 广播左矩阵第3行第3列元素
		//__m128 result3 = _mm_add_ps(
		//	_mm_add_ps(_mm_mul_ps(a3_0, mx_row0), _mm_mul_ps(a3_1, mx_row1)),
		//	_mm_add_ps(_mm_mul_ps(a3_2, mx_row2), _mm_mul_ps(a3_3, mx_row3))
		//);
		//_mm_store_ps(&ret.a[3 * 4], result3);     // 存储结果矩阵第3行

		//return ret;

		
		// -- SIMD 优化版 -- //
		matrix ret;

		__m128 mx_row0 = _mm_load_ps(&mx.a[0]);   // mx[0][0], mx[0][1], mx[0][2], mx[0][3]
		__m128 mx_row1 = _mm_load_ps(&mx.a[4]);   // mx[1][0], mx[1][1], mx[1][2], mx[1][3]
		__m128 mx_row2 = _mm_load_ps(&mx.a[8]);   // mx[2][0], mx[2][1], mx[2][2], mx[2][3]
		__m128 mx_row3 = _mm_load_ps(&mx.a[12]);  // mx[3][0], mx[3][1], mx[3][2], mx[3][3] 

		__m256 mx_2row0 = _mm256_set_m128(mx_row0, mx_row0); 
		__m256 mx_2row1 = _mm256_set_m128(mx_row1, mx_row1);
		__m256 mx_2row2 = _mm256_set_m128(mx_row2, mx_row2);
		__m256 mx_2row3 = _mm256_set_m128(mx_row3, mx_row3);

		// 两行两行计算
		//__m256 a00a10 = _mm256_set_m128(_mm_set1_ps(a[4]), _mm_set1_ps(a[0]));
		//__m256 a01a11 = _mm256_set_m128(_mm_set1_ps(a[5]), _mm_set1_ps(a[1]));
		//__m256 a02a12 = _mm256_set_m128(_mm_set1_ps(a[6]), _mm_set1_ps(a[2]));
		//__m256 a03a13 = _mm256_set_m128(_mm_set1_ps(a[7]), _mm_set1_ps(a[3]));

		/*__m256 res_row0row1 = _mm256_add_ps(
			_mm256_add_ps(_mm256_mul_ps(a00a10, mx_2row0), _mm256_mul_ps(a01a11, mx_2row1)),
			_mm256_add_ps(_mm256_mul_ps(a02a12, mx_2row2), _mm256_mul_ps(a03a13, mx_2row3))
		);*/
		__m256 res_row0row1 =
			_mm256_fmadd_ps(_mm256_set_m128(_mm_set1_ps(a[4]), _mm_set1_ps(a[0])), mx_2row0,          // a00a10*mx_2row0 + 下一个fmadd的结果
				_mm256_fmadd_ps(_mm256_set_m128(_mm_set1_ps(a[5]), _mm_set1_ps(a[1])), mx_2row1,          // a01a11*mx_2row1 + 下一个fmadd的结果
					_mm256_fmadd_ps(_mm256_set_m128(_mm_set1_ps(a[6]), _mm_set1_ps(a[2])), mx_2row2,          // a02a12*mx_2row2 + a03a13*mx_2row3
						_mm256_mul_ps(_mm256_set_m128(_mm_set1_ps(a[7]), _mm_set1_ps(a[3])), mx_2row3))));
		_mm256_storeu_ps(&ret.a[0], res_row0row1); // 存储结果的前两行

		// 第二组两行
		/*__m256 a20a30 = _mm256_set_m128(_mm_set1_ps(a[12]), _mm_set1_ps(a[8]));
		__m256 a21a31 = _mm256_set_m128(_mm_set1_ps(a[13]), _mm_set1_ps(a[9]));
		__m256 a22a32 = _mm256_set_m128(_mm_set1_ps(a[14]), _mm_set1_ps(a[10]));
		__m256 a23a33 = _mm256_set_m128(_mm_set1_ps(a[15]), _mm_set1_ps(a[11]));*/

		/*__m256 res_row2row3 = _mm256_add_ps(
			_mm256_add_ps(_mm256_mul_ps(a20a30, mx_2row0), _mm256_mul_ps(a21a31, mx_2row1)),
			_mm256_add_ps(_mm256_mul_ps(a22a32, mx_2row2), _mm256_mul_ps(a23a33, mx_2row3))
		);*/
		__m256 res_row2row3 =
			_mm256_fmadd_ps(_mm256_set_m128(_mm_set1_ps(a[12]), _mm_set1_ps(a[8])), mx_2row0,
				_mm256_fmadd_ps(_mm256_set_m128(_mm_set1_ps(a[13]), _mm_set1_ps(a[9])), mx_2row1,
					_mm256_fmadd_ps(_mm256_set_m128(_mm_set1_ps(a[14]), _mm_set1_ps(a[10])), mx_2row2,
						_mm256_mul_ps(_mm256_set_m128(_mm_set1_ps(a[15]), _mm_set1_ps(a[11])), mx_2row3))));
		_mm256_storeu_ps(&ret.a[8], res_row2row3); // 存储结果的后两行
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
		//matrix m;
		//// 展开循环
		/*m.a[0] = m.a[5] = m.a[10] = m.a[15] = 1.0f;
		m.a[1] = m.a[2] = m.a[3] = m.a[4] =
			m.a[6] = m.a[7] = m.a[8] = m.a[9] =
			m.a[11] = m.a[12] = m.a[13] = m.a[14] = 0.0f;
		return m;*/

		//-- SIMD --//
		matrix m;
		// 使用AVX指令快速设置单位矩阵
		__m256 zero = _mm256_setzero_ps();
		// __m128 one = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);  // 只有最后一个元素是1

		// 存储零向量
		_mm256_store_ps(&m.a[0], zero);
		_mm256_store_ps(&m.a[8], zero);

		// 设置对角线元素
		m.a[0] = m.a[5] = m.a[10] = m.a[15] = 1.0f;
		return m;
    }

private:
    // Set all elements of the matrix to 0
    void zero() {
		/*for (unsigned int i = 0; i < 16; i++)
			a[i] = 0.f;*/

		// unrolling
		/*a[0] = a[1] = a[2] = a[3] =
			a[4] = a[5] = a[6] = a[7] =
			a[8] = a[9] = a[10] = a[11] =
			a[12] = a[13] = a[14] = a[15] = 0.0f;*/

		//-- SIMD --//
		__m256 zero_vec = _mm256_setzero_ps();
		_mm256_store_ps(&a[0], zero_vec);
		_mm256_store_ps(&a[8], zero_vec);
    }

    // Set the matrix as an identity matrix
    void identity() {
		/*for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				m[i][j] = (i == j) ? 1.0f : 0.0f;
			}
		}*/

		// unrolling
		// 展开初始化，避免循环
		//a[0] = a[5] = a[10] = a[15] = 1.0f;
		//a[1] = a[2] = a[3] = a[4] =
		//	a[6] = a[7] = a[8] = a[9] =
		//	a[11] = a[12] = a[13] = a[14] = 0.0f;
		zero();
		a[0] = a[5] = a[10] = a[15] = 1.0f;
    }
};



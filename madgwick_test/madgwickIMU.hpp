//=====================================================================================================
// MadgwickIMU.hpp
//=====================================================================================================
//
// Implementation of Madgwick's IMU algorithm.
// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
//
// Date			Author          Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
// 19/02/2012	SOH Madgwick	Magnetometer measurement is normalised
// 19/06/2018   Juan Gallostra  Port code to c++ and update MadgwickAHRSupdateIMU signature
// 29/12/2023   Cameron King    Refactored further to c++ structures.
//=====================================================================================================

#include <cmath>

// Quaternion type

struct quat
{
	quat() {}
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
	float s = 0.0f;
	quat(float _x, float _y, float _z, float _s) : x(_x), y(_y), z(_z), s(_s) {}
};

// 3D vertex type

struct coords3d
{
	coords3d() {}
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
	coords3d(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

class madgwick
{
public:
	madgwick();
	float invSqrt(float x);
	quat MadgwickAHRSupdateIMU(coords3d gyro, coords3d accel);
	~madgwick();

private:
	float sampleFreq = 0.015f;
	volatile float beta = 0.1f;									// (betaDef) 2 * proportional gain (Kp)
	volatile float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;	// quaternion of sensor frame relative to auxiliary frame
};

madgwick::madgwick()
{

}

// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root

float madgwick::invSqrt(float x)
{
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * (1.5f - (halfx * y * y));
	return y;
}

// IMU algorithm update

quat madgwick::MadgwickAHRSupdateIMU(coords3d gyro, coords3d accel)
{
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-q1 * gyro.x - q2 * gyro.y - q3 * gyro.z);
	qDot2 = 0.5f * (q0 * gyro.x + q2 * gyro.z - q3 * gyro.y);
	qDot3 = 0.5f * (q0 * gyro.y - q1 * gyro.z + q3 * gyro.x);
	qDot4 = 0.5f * (q0 * gyro.z + q1 * gyro.y - q2 * gyro.x);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((accel.x == 0.0f) && (accel.y == 0.0f) && (accel.z == 0.0f)))
	{
		// Normalise accelerometer measurement
		recipNorm = invSqrt(accel.x * accel.x + accel.y * accel.y + accel.z * accel.z);
		accel.x *= recipNorm;
		accel.y *= recipNorm;
		accel.z *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0 = 2.0f * q0;
		_2q1 = 2.0f * q1;
		_2q2 = 2.0f * q2;
		_2q3 = 2.0f * q3;
		_4q0 = 4.0f * q0;
		_4q1 = 4.0f * q1;
		_4q2 = 4.0f * q2;
		_8q1 = 8.0f * q1;
		_8q2 = 8.0f * q2;
		q0q0 = q0 * q0;
		q1q1 = q1 * q1;
		q2q2 = q2 * q2;
		q3q3 = q3 * q3;

		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * accel.x + _4q0 * q1q1 - _2q1 * accel.y;
		s1 = _4q1 * q3q3 - _2q3 * accel.x + 4.0f * q0q0 * q1 - _2q0 * accel.y - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * accel.z;
		s2 = 4.0f * q0q0 * q2 + _2q0 * accel.x + _4q2 * q3q3 - _2q3 * accel.y - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * accel.z;
		s3 = 4.0f * q1q1 * q3 - _2q1 * accel.x + 4.0f * q2q2 * q3 - _2q2 * accel.y;
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= beta * s0;
		qDot2 -= beta * s1;
		qDot3 -= beta * s2;
		qDot4 -= beta * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	q0 += qDot1 * (1.0f / sampleFreq);
	q1 += qDot2 * (1.0f / sampleFreq);
	q2 += qDot3 * (1.0f / sampleFreq);
	q3 += qDot4 * (1.0f / sampleFreq);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;

	// store quaternion
	quat q;
	q.x = q0;
	q.y = q1;
	q.z = q2;
	q.s = q3;

	return q;
}

madgwick::~madgwick()
{

}

//====================================================================================================
// END OF CODE
//====================================================================================================
#ifndef __VECTOR_3D__
#define __VECTOR_3D__

// This file contains the defintion of the class Vector3D

#include "Matrix.h"

class Normal;
class Point3D;

//----------------------------------------- class Vector3D

class Vector3D {
	public:
	
		double	x, y, z;

	public:
	
		Vector3D(void);											// default constructor
		Vector3D(double a);										// constructor
		Vector3D(double _x, double _y, double _z);				// constructor
		Vector3D(const Vector3D& v);							// copy constructor
		Vector3D(const Normal& n);								// constructs a vector from a Normal
		Vector3D(const Point3D& p);								// constructs a vector from a point

		~Vector3D (void);										// destructor

		Vector3D& 												// assignment operator
		operator= (const Vector3D& rhs);
		
		Vector3D& 												// assign a Normal to a vector
		operator= (const Normal& rhs); 
		
		Vector3D& 												// assign a Point3D to a vector
		operator= (const Point3D& rhs); 
		
		Vector3D												// unary minus
		operator- (void) const;									
				
		double													// length
		length(void);
		
		double													// square of the length
		len_squared(void);
		
		Vector3D												// multiplication by a double on the right
		operator* (const double a) const;
				
		Vector3D												// division by a double
		operator/ (const double a) const;
		
		Vector3D												// addition
		operator+ (const Vector3D& v) const;
		
		Vector3D& 												// compound addition
		operator+= (const Vector3D& v);					
		
		Vector3D												// subtraction
		operator- (const Vector3D& v) const;					
		
		double 													// dot product							
		operator* (const Vector3D& b) const;
		
		Vector3D 												// cross product				
		operator^ (const Vector3D& v) const;
		
		void 													// convert vector to a unit vector
		normalize(void); 
		
		Vector3D& 												// return a unit vector, and normalize the vector												
		hat(void);
};


inline Vector3D 
Vector3D::operator- (void) const {
	return (Vector3D(-x, -y, -z));    
}

inline double													
Vector3D::len_squared(void) {	
	return (x * x + y * y + z * z);
}

inline Vector3D
Vector3D::operator* (const double a) const {	
	return (Vector3D(x * a, y * a, z * a));	
}

inline Vector3D
Vector3D::operator/ (const double a) const {	
	return (Vector3D(x / a, y / a, z / a));	
}


inline Vector3D
Vector3D::operator+ (const Vector3D& v) const {
	return (Vector3D(x + v.x, y + v.y, z + v.z));
}

inline Vector3D
Vector3D::operator- (const Vector3D& v) const {
	return (Vector3D(x - v.x, y - v.y, z - v.z));
}

inline double 
Vector3D::operator* (const Vector3D& v) const {
	return (x * v.x + y * v.y + z * v.z);
} 

inline Vector3D 									
Vector3D::operator^ (const Vector3D& v) const {
	return (Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x));
}

inline Vector3D& 
Vector3D::operator+= (const Vector3D& v) {
	x += v.x; y += v.y; z += v.z;
	return (*this);
}


Vector3D 
operator* (const double a, const Vector3D& v);

inline Vector3D 
operator* (const double a, const Vector3D& v) {
	return (Vector3D(a * v.x, a * v.y, a * v.z));	
}

Vector3D 
operator* (const Matrix& mat, const Vector3D& v);


#endif



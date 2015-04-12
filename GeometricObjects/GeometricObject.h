#ifndef __GEOMETRIC_OBJECT__
#define __GEOMETRIC_OBJECT__

#include <math.h>
#include "Constants.h"
#include "BBox.h"
#include "RGBColor.h"
#include "Point3D.h"
#include "Vector3D.h"
#include "Normal.h"
#include "Ray.h"
#include "ShadeRec.h"

class Material;

class GeometricObject {	
	public:	

		GeometricObject(void);
		
		GeometricObject(const GeometricObject& object);
	
		virtual GeometricObject*
		clone(void) const = 0;

		virtual
		~GeometricObject(void);
			
		virtual bool 												 
		hit(const Ray& ray, double& t, ShadeRec& s) const = 0;

		virtual void set_material(Material* mPtr);
		Material* get_material(void) const;

		void set_color(const RGBColor& c);
				
		void set_color(const float r, const float g, const float b);

		RGBColor get_color(void);

		virtual void set_bounding_box(void);

		virtual BBox get_bounding_box(void);

		virtual void add_object(GeometricObject* object_ptr);


		virtual Point3D sample(void);

		virtual float pdf(ShadeRec& sr);

		virtual Normal get_normal(void) const;

		virtual Normal get_normal(const Point3D& p);

		virtual bool
		shadow_hit(const Ray& ray, float& tmin) const;


	protected:
	
		mutable Material*   material_ptr;   	// mutable allows the const functions Compound::hit, Instance::hit, and RegularGrid::hit to assign to material_ptr
		RGBColor  color;				// only used for Bare Bones ray tracing
	
		GeometricObject& operator= (const GeometricObject& rhs);
};


inline void
GeometricObject::set_color(const RGBColor& c) {
	color = c;
}

inline void
GeometricObject::set_color(const float r, const float g, const float b) {
	color.r = r;
	color.b = b;
	color.g = g;
}


inline RGBColor GeometricObject::get_color(void) {
	return (color);
}

#endif

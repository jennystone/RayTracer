#include "PointLight.h"

// ---------------------------------------------------------------------- default constructor

PointLight::PointLight(void)
	: 	Light(),
		ls(1.0),
		color(1.0),
		location(0, 1, 0)
{}


// ---------------------------------------------------------------------- dopy constructor

PointLight::PointLight(const PointLight& dl)
	: 	Light(dl),
		ls(dl.ls),
		color(dl.color),
		location(dl.location)
{}


// ---------------------------------------------------------------------- clone

Light* 
PointLight::clone(void) const {
	return (new PointLight(*this));
}


// ---------------------------------------------------------------------- assignment operator

PointLight&
PointLight::operator= (const PointLight& rhs)
{
	if (this == &rhs)
		return (*this);
			
	Light::operator= (rhs);
	
	ls		= rhs.ls;
	color 	= rhs.color;
	location 	= rhs.location;

	return (*this);
}


// ---------------------------------------------------------------------- destructor																			

PointLight::~PointLight(void) {}


// ---------------------------------------------------------------------- get_locationection
// as this function is virtual, it shouldn't be inlined 

Vector3D								
PointLight::get_direction(ShadeRec& sr) {
	return ((location - sr.hit_point).hat());
}	

// ------------------------------------------------------------------------------  L

RGBColor
PointLight::L(ShadeRec& s) {
	return (ls * color);
}

bool
PointLight::in_shadow(const Ray& r, const ShadeRec& sr) {
	float t;
	int numObjects = sr.w.objects.size();
	float d = (location-r.o).length();

	for (int j = 0; j < numObjects; j++)
		if (sr.w.objects[j]->shadow_hit(r, t) && t < d)
			return (true);

	return (false);
}

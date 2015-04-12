#include "Directional.h"
#include <stdio.h>

// ---------------------------------------------------------------------- default constructor

Directional::Directional(void)
	: 	Light(),
		ls(1.0),
		color(1.0),
		dir(0, 1, 0)			
{}


// ---------------------------------------------------------------------- dopy constructor

Directional::Directional(const Directional& dl)
	: 	Light(dl),
		ls(dl.ls),
		color(dl.color),
		dir(dl.dir)  		
{}


// ---------------------------------------------------------------------- clone

Light* 
Directional::clone(void) const {
	return (new Directional(*this));
}


// ---------------------------------------------------------------------- assignment operator

Directional& 
Directional::operator= (const Directional& rhs) 	
{
	if (this == &rhs)
		return (*this);
			
	Light::operator= (rhs);
	
	ls		= rhs.ls;
	color 	= rhs.color;
	dir 	= rhs.dir;

	return (*this);
}


// ---------------------------------------------------------------------- destructor																			

Directional::~Directional(void) {}


// ---------------------------------------------------------------------- get_direction
// as this function is virtual, it shouldn't be inlined 

Vector3D								
Directional::get_direction(ShadeRec& sr) {
	return (dir);
}	

// ------------------------------------------------------------------------------  L

RGBColor
Directional::L(ShadeRec& s) {	
	return (ls * color);
}

bool
Directional::in_shadow(const Ray& r, const ShadeRec& sr) {
        float t;
        int num_objects = sr.w.objects.size();
        //float d = location.distance(ray.o);

        for (int j = 0; j < num_objects; j++)
                {if (sr.w.objects[j]->shadow_hit(r, t))// && t < d)
                        return (true);
                }
        return (false);
}


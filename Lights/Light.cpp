#include "Light.h"

#include "Constants.h"

// ---------------------------------------------------------------------- default constructor

Light::Light(void):
	shadow(false) {}

// ---------------------------------------------------------------------- dopy constructor

Light::Light(const Light& ls):
		shadow(ls.shadow){}


bool
Light::in_shadow(const Ray&r, const ShadeRec& sr)
{
	return false;
}
// ---------------------------------------------------------------------- assignment operator

Light& 
Light::operator= (const Light& rhs) {
	if (this == &rhs)
		return (*this);

	return (*this);
}


// ---------------------------------------------------------------------- destructor

Light::~Light(void) {} 



// ---------------------------------------------------------------------- L
// returns the radiance

RGBColor								
Light::L(ShadeRec& s) {
	return (black);
}

void
Light::set_shadow(bool b){
	shadow = b;
}

bool
Light::casts_shadow(){
	return shadow;
}

#ifndef __LIGHT__
#define __LIGHT__

#include "Vector3D.h"
#include "RGBColor.h"
#include <stdio.h>
class Ray;

class ShadeRec;


//-------------------------------------------------------------------- class Light

class Light {	
	public:
	
		Light(void);
								
		Light(const Light& ls);			

		Light& 								
		operator= (const Light& rhs); 

		virtual Light* 						
		clone(void) const = 0;
		
		virtual 							
		~Light(void);
						
		virtual bool
		in_shadow(const Ray& r,const ShadeRec& sr);

		virtual Vector3D								
		get_direction(ShadeRec& sr) = 0;				
																
		virtual RGBColor														
		L(ShadeRec& sr);

		virtual void
		set_shadow(bool b);

		virtual bool
		casts_shadow();

private:
		bool		shadow;

};

#endif

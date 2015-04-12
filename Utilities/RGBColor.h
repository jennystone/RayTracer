#ifndef __RGB_COLOR__
#define __RGB_COLOR__


class RGBColor {
	
	public:
	
		float	r, g, b;									
				
	public:
	
		RGBColor(void);										// default constructor
		RGBColor(float c);									// constructor
		RGBColor(float _r, float _g, float _b);				// constructor
		RGBColor(const RGBColor& c); 						// copy constructor
		
		~RGBColor(void);									// destructor
		
		RGBColor& 											// assignment operator
		operator= (const RGBColor& rhs); 

		RGBColor 											// addition
		operator+ (const RGBColor& c) const;	

		RGBColor& 											// compound addition
		operator+= (const RGBColor& c);
		
		RGBColor 											// multiplication by a float on the right
		operator* (const float a) const;
		
		RGBColor& 											// compound multiplication by a float on the right
		operator*= (const float a);					
				
		RGBColor 											// division by a float
		operator/ (const float a) const;
		
		RGBColor& 											// compound division by a float
		operator/= (const float a); 
				
		RGBColor 											// component-wise multiplication
		operator* (const RGBColor& c) const;
		
		bool												// are two RGBColours the same?
		operator== (const RGBColor& c) const;				

		float												// the average of the components
		average(void) const;
};


inline RGBColor 
RGBColor::operator+ (const RGBColor& c) const {
	return (RGBColor(r + c.r, g + c.g, b + c.b));
}


inline RGBColor& 
RGBColor::operator+= (const RGBColor& c) {
	r += c.r; g += c.g; b += c.b;
    return (*this);
}


inline RGBColor 
RGBColor::operator* (const float a) const {
	return (RGBColor (r * a, g * a, b * a));	
}


inline RGBColor& 
RGBColor::operator*= (const float a) {
	r *= a; g *= a; b *= a;
	return (*this);
}

inline RGBColor 
RGBColor::operator/ (const float a) const {
	return (RGBColor (r / a, g / a, b / a));
}

inline RGBColor& 
RGBColor::operator/= (const float a) {	
	r /= a; g /= a; b /= a;
	return (*this);
}

inline RGBColor 
RGBColor::operator* (const RGBColor& c) const {
	return (RGBColor (r * c.r, g * c.g, b * c.b));
} 

inline bool
RGBColor::operator== (const RGBColor& c) const {
	return (r == c.r && g == c.g && b == c.b);
}

inline float											
RGBColor::average(void) const {
	return (0.333333333333 * (r + g + b));
}


RGBColor 
operator* (const float a, const RGBColor& c);

inline RGBColor 
operator* (const float a, const RGBColor& c) {
	return (RGBColor (a * c.r, a * c.g, a * c.b));	
}


#endif

#ifndef __SHADE_REC__
#define __SHADE_REC__

// this file contains the declaration of the class ShadeRec

class World;    // only need a forward class declaration as the World data member is a reference

#include "Point3D.h"
#include "Normal.h"
#include "RGBColor.h"
#include "Ray.h"

class Material;

class ShadeRec {
        public:

                bool                            hit_an_object;          // did the ray hit an object?
                Material*                       material_ptr;           // nearest object's material
                Point3D                         hit_point;                      // world coordinates of hit point
                Point3D                         local_hit_point;        // world coordinates of hit point
                Normal                          normal;                         // normal at hit point
                RGBColor                      color;                          // used in Chapter 3 only
                Ray                                     ray;                            // for specular highlights
                float                           t;                                      // ray parameter
                int                                     depth;                          // recursion depth
                Vector3D                        dir;                            // for area lights
                World&                          w;                                      // world reference for shading
                float                                   u, v;                           // texture coordinates

                ShadeRec(World& wr);                            // constructor

                ShadeRec(const ShadeRec& sr);           // copy constructor
};

#endif

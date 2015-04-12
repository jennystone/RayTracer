#ifndef __BOX__
#define __BOX__


#include "GeometricObject.h"
#include "BBox.h"

class Box: public GeometricObject {
public:
        Box(void);

        Box(const double x0, const double x1,
                 const double y0, const double y1,
                 const double z0, const double z1);

        Box(const Point3D p0, const Point3D p1);

        Box(const Box& bbox);

        Box&
        operator=(const Box& bbox);

        virtual Box*
        clone(void) const;

        virtual bool
        hit(const Ray& ray, double& tmin, ShadeRec& sr) const;

        virtual bool
        shadow_hit(const Ray& ray, float& tmin) const;

        virtual Normal
        get_normal(const int face_hit) const;

        virtual BBox
        get_bounding_box(void) const;

private:

        double x0, x1, y0, y1, z0, z1;

};

#endif



#include "Box.h"
Box::Box (void)
: x0(-1), x1(1), y0(-1), y1(1), z0(-1), z1(1)
{}


Box::Box (const double _x0, const double _x1,
                        const double _y0, const double _y1,
                        const double _z0, const double _z1)
:       x0(_x0),
x1(_x1),
y0(_y0),
y1(_y1),
z0(_z0),
z1(_z1)
{}

Box::Box (const Point3D p0, const Point3D p1)
:       x0(p0.x),
x1(p1.x),
y0(p0.y),
y1(p1.y),
z0(p0.z),
z1(p1.z)
{}

Box::Box(const Box& bbox)
:       x0(bbox.x0),
x1(bbox.x1),
y0(bbox.y0),
y1(bbox.y1),
z0(bbox.z0),
z1(bbox.z1)
{}

Box&
Box::operator=(const Box& rhs) {
        if (this == &rhs)
                return (*this);

        GeometricObject::operator= (rhs);

        x0      = rhs.x0;
        x1      = rhs.x1;
        y0      = rhs.y0;
        y1      = rhs.y1;
        z0      = rhs.z0;
        z1      = rhs.z1;

        return (*this);
}

Box*
Box::clone(void) const {
        return (new Box(*this));
}

Normal
Box::get_normal(const int face_hit) const {
        switch(face_hit) {
                case 0:         return (Normal(-1, 0, 0));      // -x face
                case 1:         return (Normal(0, -1, 0));      // -y face
                case 2:         return (Normal(0, 0, -1));      // -z face
                case 3:         return (Normal(1, 0, 0));       // +x face
                case 4:         return (Normal(0, 1, 0));       // +y face
                case 5:         return (Normal(0, 0, 1));       // +z face
        }
}

BBox
Box::get_bounding_box(void) const {
        return (BBox(x0,x1,y0,y1,z0,z1));
}

bool
Box::hit(const Ray& ray, double& tmin, ShadeRec& sr) const {

        double ox = ray.o.x;
        double oy = ray.o.y;
        double oz = ray.o.z;

        double dx = ray.d.x;
        double dy = ray.d.y;
        double dz = ray.d.z;

        double tx_min, ty_min, tz_min;
        double tx_max, ty_max, tz_max;

        double a = 1.0 / dx;
        if (a >= 0) {
                tx_min = (x0 - ox) * a;
                tx_max = (x1 - ox) * a;
        }
        else {
                tx_min = (x1 - ox) * a;
                tx_max = (x0 - ox) * a;
        }

        double b = 1.0 / dy;
        if (b >= 0) {
                ty_min = (y0 - oy) * b;
                ty_max = (y1 - oy) * b;
        }
        else {
                ty_min = (y1 - oy) * b;
                ty_max = (y0 - oy) * b;
        }

        double c = 1.0 / dz;
        if (c >= 0) {
                tz_min = (z0 - oz) * c;
                tz_max = (z1 - oz) * c;
        }
        else {
                tz_min = (z1 - oz) * c;
                tz_max = (z0 - oz) * c;
        }

        double t0, t1;

        int face_in, face_out;

        // find largest entering t value

        if (tx_min > ty_min) {
                t0 = tx_min;
                face_in = (a >= 0.0) ? 0 : 3;
        }
        else {
                t0 = ty_min;
                face_in = (b >= 0.0) ? 1 : 4;
        }

        if (tz_min > t0) {
                t0 = tz_min;
                face_in = (c >= 0.0) ? 2 : 5;
        }

        // find smallest exiting t value

        if (tx_max < ty_max) {
                t1 = tx_max;
                face_out = (a >= 0.0) ? 3 : 0;
        }
        else {
                t1 = ty_max;
                face_out = (b >= 0.0) ? 4 : 1;
        }

        if (tz_max < t1) {
                t1 = tz_max;
                face_out = (c >= 0.0) ? 5 : 2;
        }

        if (t0 < t1 && t1 > .001) {
                if (t0 > .001) {
                        tmin = t0;                      // ray hits outside surface
                        sr.normal = get_normal(face_in);
                }
                else {
                        tmin = t1;                              // ray hits inside surface
                        sr.normal = get_normal(face_out);
                }

                sr.local_hit_point = ray.o + tmin * ray.d;
                return (true);
        }
        else
                return (false);

}



bool
Box::shadow_hit(const Ray& ray, float& tmin) const {


        double ox = ray.o.x;
        double oy = ray.o.y;
        double oz = ray.o.z;

        double dx = ray.d.x;
        double dy = ray.d.y;
        double dz = ray.d.z;

        double tx_min, ty_min, tz_min;
        double tx_max, ty_max, tz_max;

        double a = 1.0 / dx;
        if (a >= 0) {
                tx_min = (x0 - ox) * a;
                tx_max = (x1 - ox) * a;
        }
        else {
                tx_min = (x1 - ox) * a;
                tx_max = (x0 - ox) * a;
        }

        double b = 1.0 / dy;
        if (b >= 0) {
                ty_min = (y0 - oy) * b;
                ty_max = (y1 - oy) * b;
        }
        else {
                ty_min = (y1 - oy) * b;
                ty_max = (y0 - oy) * b;
        }

        double c = 1.0 / dz;
        if (c >= 0) {
                tz_min = (z0 - oz) * c;
                tz_max = (z1 - oz) * c;
        }
        else {
                tz_min = (z1 - oz) * c;
                tz_max = (z0 - oz) * c;
        }

        double t0, t1;

        int face_in, face_out;

        // find largest entering t value

        if (tx_min > ty_min) {
                t0 = tx_min;
                face_in = (a >= 0.0) ? 0 : 3;
        }
        else {
                t0 = ty_min;
                face_in = (b >= 0.0) ? 1 : 4;
        }

        if (tz_min > t0) {
                t0 = tz_min;
                face_in = (c >= 0.0) ? 2 : 5;
        }

        // find smallest exiting t value

        if (tx_max < ty_max) {
                t1 = tx_max;
                face_out = (a >= 0.0) ? 3 : 0;
        }
        else {
                t1 = ty_max;
                face_out = (b >= 0.0) ? 4 : 1;
        }

        if (tz_max < t1) {
                t1 = tz_max;
                face_out = (c >= 0.0) ? 5 : 2;
        }

        if (t0 < t1 && t1 > .001) {  // condition for a hit
                if (t0 > .001) {
                        tmin = t0;                      // ray hits outside surface
                        //sr.normal = get_normal(face_in);
                }
                else {
                        tmin = t1;                              // ray hits inside surface
                        //sr.normal = get_normal(face_out);
                }

                //sr.local_hit_point = ray.o + tmin * ray.d;
                return (true);
        }
        else
                return (false);
}

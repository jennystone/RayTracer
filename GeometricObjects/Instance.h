#ifndef __INSTANCE__
#define __INSTANCE__

/*
 *  Instance.h
 *  Ray Tracer
 *
 *  Created by NoEvilPeople on 10/15/08.
 *  Copyright 2008 jmc2385@rit.edu. All rights reserved.
 *
 */

#include "GeometricObject.h"
//#include "Box.h"

class Instance: public GeometricObject {
public:

        Instance(void);

        Instance(const GeometricObject* obj_ptr);

        Instance(const Instance& instance);

        Instance&
        operator=(const Instance& instance);

        virtual Instance*
        clone(void) const;

        void
        compute_bounding_box(void);

        virtual BBox
        get_bounding_box(void) const;

        virtual bool
        shadow_hit(const Ray& ray, float& tmin) const;

        virtual bool
        hit(const Ray& ray, double& tmin, ShadeRec& sr) const;

        void
        translate(const float dx, const float dy, const float dz);

        void
        rotate_x(const float theta);

        void
        rotate_y(const float theta);

        void
        rotate_z(const float theta);

        void
        shear(const float xy, const float xz, const float yx, const float yz, const float zx, const float zy);

        void
        scale(const float x_scale, const float y_scale, const float z_scale);

        void
        scale(const float s);

        void
        reflect_across_x_axis();

        void
        reflect_across_y_axis();

        void
        reflect_across_z_axis();

        void
        transform_texture(const bool trans);

private:

        const GeometricObject* object_ptr;      // object we're transforming
        Matrix inv_matrix;                                      // inverse of the matrix we're transforming the object with
        bool transform_the_texture;                     // whether or not to transform the texture

        static Matrix forward_matrix;           // transformation matrix
        BBox bbox;                                                      //bounding box
        //Box* box;
};

inline BBox
Instance::get_bounding_box(void) const {
        return bbox;
}

inline void
Instance::transform_texture(const bool trans) {
        transform_the_texture = trans;
}

inline void
Instance::scale(const float s) {
        scale(s,s,s);
}

#endif

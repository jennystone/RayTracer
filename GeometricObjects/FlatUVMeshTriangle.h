#ifndef __FLAT_UV_MESH_TRIANGLE__
#define __FLAT_UV_MESH_TRIANGLE__

/*
 *  FlatUVMeshTriangle.h
 *  Ray Tracer
 *
 *  Created by NoEvilPeople on 10/29/08.
 *  Copyright 2008 jmc2385@rit.edu. All rights reserved.
 *
 */


#include "FlatMeshTriangle.h"

class FlatUVMeshTriangle: public FlatMeshTriangle {
public:

        FlatUVMeshTriangle(void);

        FlatUVMeshTriangle(Mesh* _meshPtr, const int i0, const int i1, const int i2);

        virtual FlatUVMeshTriangle*
        clone(void) const;

        FlatUVMeshTriangle(const FlatUVMeshTriangle& fmt);

        virtual
        ~FlatUVMeshTriangle(void);

        FlatUVMeshTriangle&
        operator= (const FlatUVMeshTriangle& rhs);

        virtual bool
        hit(const Ray& ray, double& tmin, ShadeRec& sr) const;

        //      virtual bool
        //      shadow_hit(const Ray& ray, double& tmin) const;

};

#endif

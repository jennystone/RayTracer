

#include "Instance.h"
#include "Maths.h"
//#include "Phong.h"

Matrix
Instance::forward_matrix;       // initialize to identity matrix

Instance::Instance(void)
:       GeometricObject(),
        object_ptr(NULL),
        inv_matrix(),
        transform_the_texture(false),
        bbox(-1, 1, -1, 1, -1, 1)
{}

Instance::Instance(const GeometricObject* obj_ptr)
:       GeometricObject(),
        object_ptr(obj_ptr),
        bbox(-1, 1, -1, 1, -1, 1)
{}

Instance::Instance(const Instance& instance)
:       GeometricObject(),
        object_ptr(instance.object_ptr),
        inv_matrix(instance.inv_matrix),
        transform_the_texture(instance.transform_the_texture),
        bbox(instance.bbox)
{}

Instance&
Instance::operator=(const Instance& instance)
{
        if (this == &instance)
                return (*this);

        GeometricObject::operator= (instance);

        object_ptr = instance.object_ptr;
        inv_matrix = instance.inv_matrix;
        transform_the_texture = instance.transform_the_texture;
        bbox = instance.bbox;

        return (*this);
}

Instance*
Instance::clone(void) const {
        return (new Instance(*this));
}

void
Instance::compute_bounding_box(void) {
        BBox b = bbox;
        int size = 8;
        Point3D p[size];

        p[0] = Point3D(b.x0, b.y0, b.z0);
        p[1] = Point3D(b.x1, b.y0, b.z0);
        p[2] = Point3D(b.x0, b.y1, b.z0);
        p[3] = Point3D(b.x1, b.y1, b.z0);
        p[4] = Point3D(b.x0, b.y0, b.z1);
        p[5] = Point3D(b.x1, b.y0, b.z1);
        p[6] = Point3D(b.x0, b.y1, b.z1);
        p[7] = Point3D(b.x1, b.y1, b.z1);

//      box = new Box(Point3D(b.x0, b.y0, b.z0), Point3D(b.x1, b.y1, b.z1));
//      box->set_material(new Phong());

        double minX = kHugeValue, minY = kHugeValue, minZ = kHugeValue;
        double maxX = -kHugeValue, maxY = -kHugeValue, maxZ = -kHugeValue;

        for (int j = 0; j < size; j++)
        {
                // order shouldn't matter here, right?
                p[j] = forward_matrix * p[j];

                minX = min(minX, p[j].x);
                minY = min(minY, p[j].y);
                minZ = min(minZ, p[j].z);

                maxX = max(maxX, p[j].x);
                maxY = max(maxY, p[j].y);
                maxZ = max(maxZ, p[j].z);
        }

        bbox = BBox(Point3D(minX, minY, minZ), Point3D(maxX, maxY, maxZ));

//      box = new Box(Point3D(minX, minY, minZ), Point3D(maxX, maxY, maxZ));
//      box->set_material(new Phong());
//      Box* bot_ptr = box;

        // reset forward_matrix for next instance to use
        forward_matrix.set_identity();
}

bool
Instance::shadow_hit(const Ray& ray, float& tmin) const {
        Ray inv_ray(ray);
        inv_ray.o = inv_matrix * inv_ray.o;
        inv_ray.d = inv_matrix * inv_ray.d;

        if (object_ptr->shadow_hit(inv_ray, tmin)) {
                //sr.normal = inv_matrix * sr.normal;
                //sr.normal.normalize();

                // if object has material, use that
                if (object_ptr->get_material())
                        material_ptr = object_ptr->get_material();

//              if (!transform_the_texture)
//                      sr.local_hit_point = ray.o + t * ray.d;

                return (true);
        }
        return (false);

}

bool
Instance::hit(const Ray& ray, double& tmin, ShadeRec& sr) const {
        Ray inv_ray(ray);
        inv_ray.o = inv_matrix * inv_ray.o;
        inv_ray.d = inv_matrix * inv_ray.d;

//      Box* box_ptr = box;
//
//      if (box->hit(ray, tmin, sr)) {
//              //sr.normal = sr.normal;
//              sr.normal.normalize();
//
//              material_ptr = box->get_material();
//
//              //sr.local_hit_point = ray.o + tmin * ray.d;
//
//              return true;
//      }

        if (object_ptr->hit(inv_ray, tmin, sr)) {
                sr.normal = inv_matrix * sr.normal;
                sr.normal.normalize();

                // if object has material, use that
                if (object_ptr->get_material())
                        material_ptr = object_ptr->get_material();

                if (!transform_the_texture)
                        sr.local_hit_point = ray.o + tmin * ray.d;

                return (true);
        }
        return (false);
}

void
Instance::translate(const float dx, const float dy, const float dz) {
        Matrix inv_translation_matrix;  //temp inverse translation matrix
        Matrix translation_matrix;      // temp translation matrix

        inv_translation_matrix.m[0][3] = -dx;
        inv_translation_matrix.m[1][3] = -dy;
        inv_translation_matrix.m[2][3] = -dz;

        // post-multiply for inverse trans
        inv_matrix = inv_matrix * inv_translation_matrix;

        translation_matrix.m[0][3] = dx;
        translation_matrix.m[1][3] = dy;
        translation_matrix.m[2][3] = dz;

        forward_matrix = translation_matrix * forward_matrix;   // pre-multiply
}

void 
Instance::rotate_x(const float theta) {

        float radians = theta * PI_ON_180;

        Matrix inv_rotation_matrix;     //temp inverse rotation matrix
        Matrix rotation_matrix; // temp rotation matrix

        inv_rotation_matrix.m[1][1] = cos(radians);
        inv_rotation_matrix.m[1][2] = sin(radians);
        inv_rotation_matrix.m[2][1] = -sin(radians);
        inv_rotation_matrix.m[2][2] = cos(radians);

        // post multiply
        inv_matrix = inv_matrix * inv_rotation_matrix;

        rotation_matrix.m[1][1] = cos(radians);
        rotation_matrix.m[1][2] = -sin(radians);
        rotation_matrix.m[2][1] = sin(radians);
        rotation_matrix.m[2][2] = cos(radians);

        // pre-multiply
        forward_matrix = rotation_matrix * forward_matrix;
};

void 
Instance::rotate_y(const float theta) {

        float radians = theta * PI_ON_180;

        Matrix inv_rotation_matrix;     //temp inverse rotation matrix
        Matrix rotation_matrix; // temp rotation matrix

        inv_rotation_matrix.m[0][0] = cos(radians);
        inv_rotation_matrix.m[0][2] = -sin(radians);
        inv_rotation_matrix.m[2][0] = sin(radians);
        inv_rotation_matrix.m[2][2] = cos(radians);

        // post multiply
        inv_matrix = inv_matrix * inv_rotation_matrix;

        rotation_matrix.m[0][0] = cos(radians);
        rotation_matrix.m[0][2] = sin(radians);
        rotation_matrix.m[2][0] = -sin(radians);
        rotation_matrix.m[2][2] = cos(radians);

        // pre-multiply
        forward_matrix = rotation_matrix * forward_matrix;
};

void 
Instance::rotate_z(const float theta) {

        float radians = theta * PI_ON_180;

        Matrix inv_rotation_matrix;     //temp inverse rotation matrix
        Matrix rotation_matrix; // temp rotation matrix

        inv_rotation_matrix.m[0][0] = cos(radians);
        inv_rotation_matrix.m[0][1] = sin(radians);
        inv_rotation_matrix.m[1][0] = -sin(radians);
        inv_rotation_matrix.m[1][1] = cos(radians);

        // post multiply
        inv_matrix = inv_matrix * inv_rotation_matrix;

        rotation_matrix.m[0][0] = cos(radians);
        rotation_matrix.m[0][1] = -sin(radians);
        rotation_matrix.m[1][0] = sin(radians);
        rotation_matrix.m[1][1] = cos(radians);

        // pre-multiply
        forward_matrix = rotation_matrix * forward_matrix;
};

void
Instance::shear(const float xy, const float xz, const float yx, const float yz, const float zx, const float zy) {

        Matrix inv_shear_matrix;  // temp inverse shear matrix
        Matrix shear_matrix;    // temp sheer matrix

        inv_shear_matrix.m[0][0] = 1 - yz * zy;
        inv_shear_matrix.m[0][1] = -yx + yz * zx;
        inv_shear_matrix.m[0][2] = -zx + yx * zy;
        inv_shear_matrix.m[1][0] = -xy + xz * zy;
        inv_shear_matrix.m[1][1] = 1 - xz * zx;
        inv_shear_matrix.m[1][2] = -zy + xy * zx;
        inv_shear_matrix.m[2][0] = -xz + xy * yz;
        inv_shear_matrix.m[2][1] = -yz + xz * yx;
        inv_shear_matrix.m[2][2] = 1 - xy * yx;

        float inv_determinant = 1 / (1 - xy * yx - xz * zx - yz * zy + xy * yz * zx + xz * yx * zy);

        // post multiply
        inv_matrix = inv_matrix * (inv_shear_matrix.scalar_mult(inv_determinant));

        shear_matrix.m[0][1] = yx;
        shear_matrix.m[0][2] = zx;
        shear_matrix.m[1][0] = xy;
        shear_matrix.m[1][2] = zy;
        shear_matrix.m[2][0] = xz;
        shear_matrix.m[2][1] = yz;

        // pre-multiply
        forward_matrix = shear_matrix * forward_matrix;
}

void
Instance::scale(const float x_scale, const float y_scale, const float z_scale) {

        Matrix inv_scale_matrix;        // temp inverse scale matrix
        Matrix scale_matrix;    // temp scale matrix

        inv_scale_matrix.m[0][0] = 1 / x_scale;
        inv_scale_matrix.m[1][1] = 1 / y_scale;
        inv_scale_matrix.m[2][2] = 1 / z_scale;

        // post-multiply
        inv_matrix = inv_matrix * inv_scale_matrix;

        scale_matrix.m[0][0] = x_scale;
        scale_matrix.m[1][1] = y_scale;
        scale_matrix.m[2][2] = z_scale;

        // pre multiply
        forward_matrix = scale_matrix * forward_matrix;
}

void
Instance::reflect_across_x_axis() {

        Matrix inv_reflect_matrix;      // temp inverse reflect matrix
        Matrix reflect_matrix;  // temp reflect matrix

        inv_reflect_matrix.m[0][0] = -1;

        // post multiply
        inv_matrix = inv_matrix * inv_reflect_matrix;

        reflect_matrix.m[1][1] = -1;
        reflect_matrix.m[2][2] = -1;

        // pre-multiply
        forward_matrix = reflect_matrix * forward_matrix;
}

void
Instance::reflect_across_y_axis() {
        Matrix inv_reflect_matrix;      // temp inverse reflect matrix
        Matrix reflect_matrix;  // temp reflect matrix

        inv_reflect_matrix.m[1][1] = -1;

        // post multiply
        inv_matrix = inv_matrix * inv_reflect_matrix;

        reflect_matrix.m[0][0] = -1;
        reflect_matrix.m[2][2] = -1;

        // pre-multiply
        forward_matrix = reflect_matrix * forward_matrix;
}

void
Instance::reflect_across_z_axis(){
        Matrix inv_reflect_matrix;      // temp inverse reflect matrix
        Matrix reflect_matrix;  // temp reflect matrix

        inv_reflect_matrix.m[2][2] = -1;

        // post multiply
        inv_matrix = inv_matrix * inv_reflect_matrix;

        reflect_matrix.m[0][0] = -1;
        reflect_matrix.m[1][1] = -1;

        // pre-multiply
        forward_matrix = reflect_matrix * forward_matrix;
}

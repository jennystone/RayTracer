#ifndef __GRID__
#define __GRID__


#include "Compound.h"
#include <vector>
#include "Mesh.h"

class Grid: public Compound {
public:

        Grid(void);

        Grid(const Grid& grid);

        Grid&
        operator=(const Grid& grid);

        Grid(Mesh* _mesh_ptr);

        ~Grid(void);

        virtual Grid*                                                           // virtual copy constructor
        clone(void) const;

        virtual BBox
        get_bounding_box(void) const;

        void
        setup_cells(void);

        virtual bool
        hit(const Ray& ray, double& t, ShadeRec& s) const;

        virtual bool
        shadow_hit(const Ray& ray, float& tmin) const;

        void
        read_flat_triangles(char* file_name);

        void
        read_smooth_triangles(char* file_name);

        void
        reverse_mesh_normals(void);

//      void
//      tessellate_flat_sphere(const int horizontal_steps, const int vertical_steps);
//
//      void
//      tessellate_smooth_sphere(const int horizontal_steps, const int vertical_steps);

//      virtual void
//      set_material(Material* mat_ptr);

        void
        read_flat_uv_triangles(char* file_name);


        // ----------------------------------------------------------------------------- read_smooth_uv_triangles

        void
        read_smooth_uv_triangles(char* file_name);

        void
        read_uv_ply_file(char* file_name, const int triangle_type);

private:

        std::vector<GeometricObject*> cells;            // 1D array
        BBox bbox;
        int nx, ny, nz;                 // number of cells in the x, y, and z directions
        Mesh*                                           mesh_ptr;               // holds triangle data
        bool                                            reverse_normal; // some PLY files have normals that point inwards

        Point3D
        min_coordinates(void);

        Point3D
        max_coordinates(void);

        void
        delete_cells(void);

        void
        copy_cells(const std::vector<GeometricObject*>& rhs_cells);

        void
        read_ply_file(char* file_name, const int triangle_type);

        void
        compute_mesh_normals(void);

};

inline void
Grid::reverse_mesh_normals(void) {
        reverse_normal = true;
}

#endif

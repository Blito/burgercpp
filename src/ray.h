#ifndef RAY_H
#define RAY_H

#include "LinearMath/btVector3.h"

class material;

namespace ray_physics {

struct ray
{
    btVector3 from, direction;
    unsigned int depth;
    const material & media;
    float intensity, frequency;
    unsigned int distance_traveled; // [mm]
    unsigned short parent_collision; // position in collision vector

    static constexpr unsigned int max_depth = 5;
    static constexpr float intensity_epsilon = 1e-8;
};

struct segment
{
    btVector3 from, to, direction;
    float initial_intensity, attenuation;
    unsigned int distance_traveled; // [mm]
};

struct collision
{
    btVector3 position;
    unsigned short parent_collision; // position in collision vector
};

struct hit_result { ray reflection, refraction; };

hit_result hit_boundary(const ray & r, const btVector3 &hit_point, const btVector3 & surface_normal, const material & media);

// Advance through homogeneous media and decrease intensity accordingly
void travel(ray & r, float mm);

bool should_travel(const ray & r);

float max_ray_length(const ray & r);

btVector3 snells_law(const btVector3 & ray_direction, const btVector3 & surface_normal, float refr_index_1, float refr_index_2);

float reflected_intensity(float intensity_in, float media_1, float incidence_angle, float media_2, float refracted_angle);

} // end ray_physics

#endif // RAY_H

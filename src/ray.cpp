#include "ray.h"

#include <cmath>

#include "material.h"

using namespace ray_physics;

ray_physics::hit_result ray_physics::hit_boundary(const ray & r, const btVector3 & hit_point, const btVector3 & surface_normal, const material & media)
{
    const auto refraction_direction = snells_law(r.direction, surface_normal, r.media.impedance, media.impedance);

    auto refraction_angle = refraction_direction.dot(-surface_normal);
    if (refraction_angle < 0)
    {
        refraction_angle = refraction_direction.dot(surface_normal);
    }

    btScalar incidence_angle = r.direction.dot(-surface_normal);
    // TODO: Check if this if is needed.
//    if (incidence_angle < 0)
//    {
//        incidence_angle = r.direction.dot(surface_normal);
//    }

    const btVector3 reflection_direction = r.direction + 2*incidence_angle * surface_normal;
    const auto intensity_refl = reflected_intensity(r.intensity,
                                                    r.media.impedance, incidence_angle,
                                                    media.impedance, refraction_angle);
    const auto intensity_refr = r.intensity - intensity_refl;

    // Add two more rays to the stack
    ray refraction_ray { hit_point, refraction_direction, r.depth+1, media, intensity_refr > ray::intensity_epsilon ? intensity_refr : 0.0f, r.frequency, r.distance_traveled, 0 };

    ray reflection_ray { hit_point, reflection_direction, r.depth+1, r.media, intensity_refl > ray::intensity_epsilon ? intensity_refl : 0.0f, r.frequency, r.distance_traveled, 0 };

    return { reflection_ray, refraction_ray };
}

void ray_physics::travel(ray & r, float mm)
{
    r.distance_traveled += mm;
    r.intensity = r.intensity * std::exp(-r.media.attenuation*(mm*0.1f)*r.frequency);
}

bool ray_physics::should_travel(const ray & r)
{
    return r.depth < r.max_depth;
}

float ray_physics::max_ray_length(const ray & r)
{
    return 10.f /*<- cm to mm*/ * std::log(ray::intensity_epsilon/r.intensity) / -r.media.attenuation * r.frequency;
}

btVector3 ray_physics::snells_law(const btVector3 & ray_direction, const btVector3 & surface_normal, float refr_index_1, float refr_index_2)
{
    // For more details, read https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
    float cos_sigma1 = ray_direction.dot(-surface_normal);
    if (cos_sigma1 < 0.0f)
    {
        cos_sigma1 = ray_direction.dot(surface_normal);
    }

    //btVector3 reflection_dir{ 1 + 2 * cos_sigma1 * surface_normal };

    const btVector3 & l = ray_direction;
    const btVector3 & n = surface_normal;
    const float c = cos_sigma1;
    const float r = refr_index_1 / refr_index_2;

    return btVector3( r * l + (r*c - std::sqrt(1 - r*r * (1 - c*c))) * n );
}

float ray_physics::reflected_intensity(float intensity_in, float media_1, float incidence_angle, float media_2, float refracted_angle)
{
    auto && num = media_1 * incidence_angle - media_2 * refracted_angle;
    auto && denom = media_1 * incidence_angle + media_2 * refracted_angle;

    return intensity_in * pow(num/denom, 2);
}

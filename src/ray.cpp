#include "ray.h"

#include <cmath>
#include <iostream>

#include "mesh.h"

using namespace ray_physics;

ray_physics::hit_result ray_physics::hit_boundary(const ray & r, const btVector3 & hit_point, const btVector3 & surface_normal, const mesh & collided_mesh)
{
    btScalar incidence_angle = r.direction.dot(-surface_normal); // cos theta_1

    // TODO: this logic can probably be simpler
    const material * material_after_vascularities = nullptr;
    const auto & material_after_collision = [&r, &collided_mesh, &material_after_vascularities]() -> const material &
    {
        if (r.media_outside) // if we are in a vessel
        {
            if (collided_mesh.is_vascular) // and we collided against a vessel (assuming same vessel, so we're getting out of it)
            {
                material_after_vascularities = nullptr;
                return *r.media_outside; // we are going back to the stored media
            }
            else // we are still inside the vessel but went out of the surrounding organ
            {
                // update the surrounding tissue
                material_after_vascularities = r.media_outside == &collided_mesh.material_inside ? &collided_mesh.material_outside : &collided_mesh.material_inside;

                // but we remain in the same media, i.e. the vessel
                return r.media;
            }
        }
        else // we are not in a vessel
        {
            if (collided_mesh.is_vascular) // and we collided with a vessel
            {
                // update the surrounding tissue
                material_after_vascularities = &r.media; // we will come back to this tissue after getting out of the vessel
                return collided_mesh.material_inside;
            }
            else // and we collided with a regular organ
            {
                material_after_vascularities = nullptr;
                return &r.media == &collided_mesh.material_inside ? collided_mesh.material_outside : collided_mesh.material_inside;
            }
        }
    }();

    if (incidence_angle < 0)
    {
        incidence_angle = r.direction.dot(surface_normal);
    }

    const float refr_ratio = r.media.impedance / material_after_collision.impedance;
    float refraction_angle = 1 - refr_ratio*refr_ratio * (1 - incidence_angle*incidence_angle);
    const bool total_internal_reflection = refraction_angle < 0;

    refraction_angle = std::sqrt(refraction_angle);

    const auto refraction_direction = snells_law(r.direction, surface_normal, incidence_angle, refraction_angle, refr_ratio);
    const btVector3 reflection_direction = r.direction + 2*incidence_angle * surface_normal;

    const auto intensity_refl = total_internal_reflection ?
                                    r.intensity :
                                    reflection_intensity(r.intensity,
                                        r.media.impedance, incidence_angle,
                                        material_after_collision.impedance, refraction_angle);
    const auto intensity_refr = r.intensity - intensity_refl;

    // Eq. 10 in Burger13
    const float back_to_transducer_intensity = reflected_intensity(r.intensity, incidence_angle, r.media, material_after_collision);

    // Add two more rays to the stack
    ray refraction_ray { hit_point, refraction_direction, r.depth+1, material_after_collision, material_after_vascularities, intensity_refr > ray::intensity_epsilon ? intensity_refr : 0.0f, r.frequency, r.distance_traveled, 0 };

    ray reflection_ray { hit_point, reflection_direction, r.depth+1, r.media, r.media_outside, intensity_refl > ray::intensity_epsilon ? intensity_refl : 0.0f, r.frequency, r.distance_traveled, 0 };

    return { back_to_transducer_intensity, reflection_ray, refraction_ray };
}

void ray_physics::travel(ray & r, units::length::millimeter_t mm)
{
    r.distance_traveled = r.distance_traveled + mm;
    r.intensity = r.intensity * std::exp(-r.media.attenuation*(mm.to<float>()*0.01f)*r.frequency); // TODO: that 0.01 should be 0.1
}

bool ray_physics::should_travel(const ray & r)
{
    return r.depth < r.max_depth;
}

float ray_physics::max_ray_length(const ray & r)
{
    return 10.f /*<- cm to mm*/ * std::log(ray::intensity_epsilon/r.intensity) / -r.media.attenuation * r.frequency;
}

btVector3 ray_physics::snells_law(const btVector3 & ray_direction, const btVector3 & surface_normal, float incidence_angle, float refraction_angle, float refr_ratio)
{
    // For more details, read https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
    const btVector3 & l = ray_direction;
    const btVector3 & n = surface_normal;
    const float c = incidence_angle;
    const float r = refr_ratio;

    return btVector3( r * l + (r*c - refraction_angle) * n );
}

float ray_physics::reflection_intensity(const float intensity_in, const float media_1, const float incidence_angle, const float media_2, const float refracted_angle)
{
    const auto && num = media_1 * incidence_angle - media_2 * refracted_angle;
    const auto && denom = media_1 * incidence_angle + media_2 * refracted_angle;

    return intensity_in * pow(num/denom, 2);
}

float ray_physics::reflected_intensity(const float ray_intensity, const float incidence_angle, const material & ray_media, const material & colliding_media)
{
    // Eq. 10 in Burger13
    constexpr auto small_reflections_enhancement_factor = 0.2;

    constexpr auto custom_reflection_enhancement_factor = 0.05; // we made this up

    const auto specular_factor = std::pow(incidence_angle, colliding_media.specularity);
    const auto impedance_factor = std::pow(( (colliding_media.impedance - ray_media.impedance)
                                            /(colliding_media.impedance + ray_media.impedance)),2);
    const auto intensity = std::pow(ray_intensity, small_reflections_enhancement_factor);

    //std::cout << media_1.impedance << " " << media_2.impedance << std::endl;
    //std::cout << ray_intensity << ", " << specular_factor << " * " << impedance_factor << " * " << intensity << " = " << std::abs(specular_factor * impedance_factor * intensity) << std::endl;
    //return std::pow(std::abs(specular_factor * impedance_factor * ray_intensity), small_reflections_enhancement_factor);
    return std::abs(specular_factor * std::pow(impedance_factor, custom_reflection_enhancement_factor) * intensity);
}

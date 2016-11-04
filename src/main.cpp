#include "scene.h"
#include "volume.h"
#include "psf.h"
#include "rfimage.h"

#include <cmath>

constexpr unsigned int speed_of_sound = 1500; // [μm/μs], [m/s]
constexpr float transducer_frequency = 4.5f; // [Mhz]
constexpr float axial_resolution = 1.45f / transducer_frequency; // [mm], the division can be deduced from Burger13
constexpr unsigned int transducer_elements = 256;
constexpr unsigned int ultrasound_depth = 150000; // [15cm -> μm]
constexpr unsigned int max_travel_time = ultrasound_depth / speed_of_sound; // [μs]

constexpr unsigned int resolution = 145; // [μm], from Burger13
using psf_ = psf<13, 7, 7, resolution>;
using volume_ = volume<256, resolution>;
using rf_image_ = rf_image<transducer_elements, max_travel_time, static_cast<unsigned int>(axial_resolution*1000.0f/*mm->μm*/)>;

namespace
{
    float convolution(const volume_ & v, const psf_ & p, const float x, const float y, const float z)
    {
        return 0.0f;
    }
}

int main(int argc, char** argv)
{

    static const volume_ volume;
    const psf_ psf(transducer_frequency, 0.1f, 0.3f, 0.4f);

    rf_image_ rf_image;

    scene scene;
    scene.init();

    scene.step(1000.0f);
    while(true)
    {
        auto rays = scene.cast_rays<transducer_elements>();

        for (unsigned int ray_i = 0; ray_i < rays.size(); ray_i++)
        {
            const auto & ray = rays[ray_i];
            for (auto & segment : ray)
            {
                const auto steps = scene.distance(segment.from, segment.to) / axial_resolution;
                const auto delta_step = axial_resolution * segment.direction;
                const auto time_step = max_travel_time / steps; // [μs]

                auto point = segment.from;
                auto time_elapsed = 0.0f;
                auto intensity = segment.initial_intensity;
                for (unsigned int step = 0; step < steps; step++)
                {
                    float echo = intensity * convolution(volume, psf, point.getX(), point.getY(), point.getZ());

                    rf_image.add_echo(ray_i, echo, time_elapsed);

                    // Step forward through the segment, decreasing intensity using Beer-Lambert's law
                    point += delta_step;
                    time_elapsed += time_step;
                    intensity *= std::exp(-segment.attenuation * axial_resolution * transducer_frequency);
                }
            }
        }
    }

}

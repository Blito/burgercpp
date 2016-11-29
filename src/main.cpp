#include "scene.h"
#include "volume.h"
#include "psf.h"
#include "rfimage.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <units/units.h>
#include <nlohmann/json.hpp>

using namespace units::literals;
using namespace units::velocity;
using namespace units::length;
using namespace units::time;

constexpr meters_per_second_t speed_of_sound = 1500_m / 1_s; // [μm/μs], [m/s]
constexpr float transducer_frequency = 4.5f; // [Mhz]
constexpr millimeter_t axial_resolution = millimeter_t(1.45f / transducer_frequency); // [mm], the division can be deduced from Burger13
constexpr unsigned int transducer_elements = 256;
constexpr centimeter_t ultrasound_depth = 15_cm; // [15cm -> μm]
constexpr microsecond_t max_travel_time = microsecond_t(ultrasound_depth / speed_of_sound); // [μs]

constexpr unsigned int resolution = 145; // [μm], from Burger13
using psf_ = psf<13, 7, 7, resolution>;
using volume_ = volume<256, resolution>;
using rf_image_ = rf_image<transducer_elements, max_travel_time.to<unsigned int>(), static_cast<unsigned int>(axial_resolution.to<float>()*1000.0f/*mm->μm*/)>;

namespace
{
    float convolution(const volume_ & v, const material & m, const psf_ & p, const float x, const float y, const float z)
    {
        float total = 0.0f;

        for (size_t i = 0; i < p.get_axial_size(); i++)
        {
            float x_volume = x + i - p.get_axial_size()/2 * v.get_resolution_in_millis();

            for (size_t j = 0; j < p.get_lateral_size(); j++)
            {
                float y_volume = y + j - p.get_lateral_size()/2 * v.get_resolution_in_millis();

                for (size_t k = 0; k < p.get_elevation_size(); k++)
                {
                    float z_volume = z + k - p.get_elevation_size()/2 * v.get_resolution_in_millis();

                    total += v.get_scattering(m.mu1, m.mu0, m.sigma, x_volume, y_volume, z_volume) * p.get(i,j,k);
                }
            }
        }

        return total;
    }
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cout << "Incorrect argument list." << std::endl;
        return 0;
    }

    static const volume_ texture_volume;

    const psf_ psf { transducer_frequency, 0.1f, 0.3f, 0.4f };

    rf_image_ rf_image;

    nlohmann::json json;
    {
        std::ifstream infile { argv[1] };

        json << infile;
    }

    try
    {
        scene scene { json };
        scene.init();

        scene.step(1000.0f);
        while(true)
        {
            rf_image.clear();

            auto rays = scene.cast_rays<transducer_elements>();

            for (unsigned int ray_i = 0; ray_i < rays.size(); ray_i++)
            {
                const auto & ray = rays[ray_i];
                for (auto & segment : ray)
                {
                    const auto starting_micros = rf_image.micros_traveled(segment.distance_traveled /*mm -> μm*/);
                    const auto distance = scene.distance(segment.from, segment.to); // [mm]
                    const auto steps = distance / axial_resolution;
                    const auto delta_step = axial_resolution.to<float>() * segment.direction;
                    const auto time_step = rf_image.micros_traveled(axial_resolution); // [μs]

                    auto point = segment.from;
                    auto time_elapsed = starting_micros;
                    auto intensity = segment.initial_intensity;

                    for (unsigned int step = 0; step < steps && time_elapsed < max_travel_time; step++)
                    {
                        float echo = intensity * convolution(texture_volume, segment.media, psf, point.getX(), point.getY(), point.getZ());

                        rf_image.add_echo(ray_i, echo, time_elapsed);

                        // Step forward through the segment, decreasing intensity using Beer-Lambert's law
                        point += delta_step;
                        time_elapsed = time_elapsed + time_step;

                        constexpr auto k = 0.05f;
                        intensity *= std::exp(-segment.attenuation * axial_resolution.to<float>()*0.1f * transducer_frequency * k);
                    }

                    // Add reflection term, i.e. intensity directly reflected back to the transducer. See Burger13, Eq. 10.
                    rf_image.add_echo(ray_i, segment.reflected_intensity, starting_micros + time_step * (steps-1));

                }

            }

            rf_image.envelope();

            rf_image.show();
        }
    }
    catch (const std::exception & ex)
    {
        std::cout << "The program found an error and will terminate.\n"
                  << "Reason:\n"
                  << ex.what() << std::endl;
    }

}

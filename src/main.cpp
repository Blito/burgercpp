#include "scene.h"
#include "volume.h"
#include "psf.h"
#include "rfimage.h"
#include "transducer.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <units/units.h>
#include <nlohmann/json.hpp>

using namespace units::literals;
using namespace units::velocity;
using namespace units::length;
using namespace units::time;
using namespace units::angle;

constexpr meters_per_second_t speed_of_sound = 1500_m / 1_s; // [μm/μs], [m/s]
constexpr float transducer_frequency = 4.5f; // [Mhz]
constexpr millimeter_t axial_resolution = millimeter_t(1.45f / transducer_frequency); // [mm], the division can be deduced from Burger13
constexpr size_t transducer_elements = 512;
constexpr radian_t transducer_amplitude = 60_deg;
constexpr centimeter_t transducer_radius = 3_cm;
constexpr centimeter_t ultrasound_depth = 15_cm; // [15cm -> μm]
constexpr microsecond_t max_travel_time = microsecond_t(ultrasound_depth / speed_of_sound); // [μs]

constexpr unsigned int resolution = 145;//145; // [μm], from Burger13
using psf_ = psf<7, 13, 7, resolution>;
using volume_ = volume<256, resolution>;
using rf_image_ = rf_image<transducer_elements, max_travel_time.to<unsigned int>(), static_cast<unsigned int>(axial_resolution.to<float>()*1000.0f/*mm->μm*/)>;
using transducer_ = transducer<transducer_elements>;

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cout << "Incorrect argument list." << std::endl;
        return 0;
    }

    static const volume_ texture_volume;

    const psf_ psf { transducer_frequency, 0.05f, 0.2f, 0.1f };

    rf_image_ rf_image { transducer_radius, transducer_amplitude };

    nlohmann::json json;
    {
        std::ifstream infile { argv[1] };

        json << infile;
    }

    const auto & t_pos = json.at("transducerPosition");
    millimeter_t transducer_element_separation = transducer_amplitude.to<float>() * transducer_radius / transducer_elements;

    const auto & t_dir = json.at("transducerAngles");
    std::array<units::angle::degree_t, 3> transducer_angles = {degree_t((float)t_dir[0]), degree_t((float)t_dir[1]), degree_t((float)t_dir[2])};

    transducer_ transducer(transducer_frequency, transducer_radius, transducer_element_separation,
                           btVector3(t_pos[0], t_pos[1], t_pos[2]), transducer_angles);
    //btVector3(-17.0, 1.2, 6.45) // liver
    //btVector3(-19.5, 1.2, -0.45) // kidney

    std::cout << max_travel_time << std::endl;

    try
    {
        scene scene { json, transducer };

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
                        float scattering = texture_volume.get_scattering(segment.media.mu1, segment.media.mu0, segment.media.sigma, point.x(), point.y(), point.z());

                        rf_image.add_echo(ray_i, intensity * scattering, time_elapsed);

                        // Step forward through the segment, decreasing intensity using Beer-Lambert's law
                        point += delta_step;
                        time_elapsed = time_elapsed + time_step;

                        constexpr auto k = 1.0f;
                        intensity *= std::exp(-segment.attenuation * axial_resolution.to<float>()*0.01f * transducer_frequency * k);
                    }

                    // Add reflection term, i.e. intensity directly reflected back to the transducer. See Burger13, Eq. 10.
                    rf_image.add_echo(ray_i, segment.reflected_intensity, starting_micros + time_step * (steps-1));

                }

            }

            rf_image.convolve(psf);

            rf_image.envelope();

            rf_image.postprocess();

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

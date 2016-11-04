#ifndef RFIMAGE_H
#define RFIMAGE_H

#include <array>

template <unsigned int columns, unsigned int max_travel_time /*μs*/, unsigned int axial_resolution /*μm*/, unsigned int speed_of_sound = 1500 /*μm/μs*/>
class rf_image
{
public:
    void add_echo(const unsigned int column, const float echo_intensity, const float micros_from_source)
    {
        const auto row = static_cast<unsigned int>(micros_from_source / (axial_resolution / speed_of_sound));
        if (row < max_rows)
        {
            image[column][row] += echo_intensity;
        }
    }

private:
    static constexpr unsigned int max_rows = (max_travel_time * speed_of_sound) / axial_resolution;

    std::array<std::array<float, max_rows>, columns> image {{{{0.0f}}}};
};

#endif // RFIMAGE_H

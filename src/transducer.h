#ifndef TRANSDUCER_H
#define TRANSDUCER_H

#include <units/units.h>
#include <LinearMath/btVector3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

#define M_PI 3.14159

template<size_t transducer_elements>
class transducer
{
public:
    struct transducer_element
    {
        btVector3 position;
        btVector3 direction;
    };

    transducer(const float frequency, const units::length::centimeter_t radius, units::length::millimeter_t transducer_element_separation,
               const btVector3 position, const btVector3 direction) :
        frequency(frequency),
        radius(radius),
        position(position),
        direction(direction)
    {
        using namespace units::angle;
        using namespace units::literals;

        assert(transducer_element_separation * transducer_elements < M_PI * radius);

        radian_t starting_angle { 90_deg };

        auto amp = transducer_element_separation / radius;
        const radian_t amplitude { amp.to<float>() }; // angle covered by a single TE
        const radian_t angle_center_of_element { amplitude / 2.0f };

        radian_t angle = starting_angle + -(amplitude * transducer_elements / 2) + angle_center_of_element;

        for (size_t t = 0; t < transducer_elements; t++)
        {
            elements[t] = transducer_element
            {
                position + radius.to<float>() * btVector3 ( std::sin(angle.to<float>()), 0, std::cos(angle.to<float>()) ), // position
                btVector3 ( std::sin(angle.to<float>()), 0, std::cos(angle.to<float>()) )  // direction
            };

            angle = angle + amplitude;
        }

    }

    transducer_element element(size_t i) const
    {
        return elements.at(i);
    }

    void print(bool direction) const
    {
        auto print_vec = [](const auto & v)
        {
            std::cout << v.x() << "," << v.z() << std::endl;
        };

        for (auto & element : elements)
        {
            print_vec(direction? element.direction : element.position);
        }
    }

    const float frequency;

    const btVector3 position, direction;

private:
    const units::length::millimeter_t radius;

    std::array<transducer_element, transducer_elements> elements;
};

#endif // TRANSDUCER_H

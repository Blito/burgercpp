#ifndef PSF_H
#define PSF_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include <ratio>

#define M_PI 3.14159

/**
 * The PSF has values in a voxelized space.
 *
 * It is defined amongst a bounded range, since outside of it it's value is 0.
 *
 * It has three diferent ranges: axial, lateral and elevation.
 * The axial range varies according to frequency.
 * Lateral and elevation ranges vary according to distance to the transducer.
 *
 * The discretization is done around the center of the psf. Each voxel has a
 * constant resolution in mm. This means that the size of the matrix that holds
 * the discrete volume should be able to hold every value of the psf while the
 * psf' ranges are at a maximum. This also means that when the psf' ranges become
 * smaller (near the focus zone), there will be zeroed voxels.
 */
template <size_t axial_size, size_t lateral_size, size_t elevation_size, unsigned int resolution_micrometers>
class psf
{
    static_assert(axial_size % 2, "axial_size must be an odd positive integer");
    static_assert(lateral_size % 2, "lateral_size must be an odd positive integer");
    static_assert(elevation_size % 2, "elevation_size must be an odd positive integer");

public:
    psf(const float freq, const float var_x, const float var_y, const float var_z) :
        freq(freq),
        var_x(var_x),
        var_y(var_y),
        var_z(var_z)
    {
        constexpr auto half_axial = axial_size * resolution_micrometers / 1000.0f / 2.0f; // [mm]
        constexpr auto half_lateral = lateral_size * resolution_micrometers / 1000.0f / 2.0f; // [mm]
        constexpr auto half_elevation = elevation_size * resolution_micrometers / 1000.0f / 2.0f; // [mm]
        constexpr auto resolution = resolution_micrometers / 1000.0f; // [mm]

        // Fill axial kernel
        for (size_t i = 0; i < axial_size; i++)
        {
            const float x = i * resolution - half_axial; // [mm]
            axial_kernel[i] = axial_function(x);
        }

        // Fill lateral kernel
        for (size_t i = 0; i < lateral_size; i++)
        {
            const float y = i * resolution - half_lateral; // [mm]
            lateral_kernel[i] = lateral_function(y);
        }
    }

    constexpr size_t get_axial_size() const
    {
        return axial_size;
    }

    constexpr size_t get_lateral_size() const
    {
        return lateral_size;
    }

    constexpr size_t get_elevation_size() const
    {
        return elevation_size;
    }

    std::array<float, axial_size> axial_kernel;
    std::array<float, lateral_size> lateral_kernel;
    std::array<float, elevation_size> elevation_kernel;

private:
    float axial_function(const float x) const
    {
        using namespace std;

        return exp(-0.5f*(pow(x,2)/var_x))*cos(2*M_PI*freq*x);
    }

    float lateral_function(const float y) const
    {
        using namespace std;

        return exp(-0.5f*(pow(y,2)/var_y));
    }

    const float var_x, var_y, var_z;
    const float freq;
};

#endif // PSF_H

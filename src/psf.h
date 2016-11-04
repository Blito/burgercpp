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

        for (size_t i = 0; i < axial_size; i++)
        {
            const float x = i * resolution - half_axial; // [mm]
            for (size_t j = 0; j < lateral_size; j++)
            {
                const float y = j * resolution - half_lateral; // [mm]
                for (size_t k = 0; k < elevation_size; k++)
                {
                    const float z = k * resolution - half_elevation; // [mm]
                    values[i][j][k] = function(x,y,z);
                }
            }
        }
    }

    float get(const unsigned int x, const unsigned int y, const unsigned int z) const
    {
        return values[x][y][z];
    }

private:
    float function(const float x, const float y, const float z) const
    {
        using namespace std;

        return exp(-0.5f*(pow(x,2)/var_x + pow(y,2)/var_y + pow(z,2)/var_z))*cos(2*M_PI*freq*x);
    }

    std::array<std::array<std::array<float, elevation_size>, lateral_size>, axial_size> values;

    const float var_x, var_y, var_z;
    const float freq;
};

#endif // PSF_H

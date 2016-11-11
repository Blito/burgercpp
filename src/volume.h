#ifndef VOLUME_H
#define VOLUME_H

#include <array>
#include <random>

/**
 * The volume class holds a cubic matrix with gaussian random values where
 * each voxel has a mm resolution.
 *
 * It is posible to query the volume for one of its values. The volume loops
 * over its internal matrix and returns the expected value.
 */
template <unsigned int size, unsigned int resolution_micrometers>
class volume
{
public:
    volume()
    {
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0,1.0);

        for (unsigned int i = 0; i < size; i++)
        {
            for (unsigned int j = 0; j < size; j++)
            {
                for (unsigned int k = 0; k < size; k++)
                {
                    matrix[i][j][k].texture_noise = distribution(generator);
                    matrix[i][j][k].scattering_probability = distribution(generator);
                }
            }
        }
    }

    constexpr float get_resolution_in_millis() const
    {
        return static_cast<float>(resolution_micrometers)/1000.0f;
    }

    /**
     * Gets the scattering value for a tissue with given properties at a fixed point in space.
     * @sa Eq. 15 in Burger13
     */
    float get_scattering(const float scattering_density, const float scattering_mu, const float scattering_sigma,
              const float x_millis, const float y_millis, const float z_millis) const
    {
        constexpr float resolution = resolution_micrometers / 1000.0f; // [mm]

        // TODO: Change static_cast to linear interpolation?
        const unsigned int x = static_cast<unsigned int>(x_millis / resolution) % size;
        const unsigned int y = static_cast<unsigned int>(y_millis / resolution) % size;
        const unsigned int z = static_cast<unsigned int>(z_millis / resolution) % size;

        const auto & voxel = matrix[x][y][z];

        return voxel.scattering_probability >= scattering_density ?
                    voxel.texture_noise * scattering_sigma + scattering_mu :
                    0.0f;
    }

private:
    struct voxel_values
    {
        float texture_noise;
        float scattering_probability;
    };

    std::array<std::array<std::array<voxel_values, size>, size>, size> matrix;
};

#endif // VOLUME_H

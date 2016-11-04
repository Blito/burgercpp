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
//        std::default_random_engine generator;
//        std::normal_distribution<double> distribution(0.0,1.0);

//        for (unsigned int i = 0; i < size; i++)
//        {
//            for (unsigned int j = 0; j < size; j++)
//            {
//                for (unsigned int k = 0; k < size; k++)
//                {
//                    matrix[i][j][k] = distribution(generator);
//                }
//            }
//        }
    }

    float get(const float x_millis, const float y_millis, const float z_millis) const
    {
        constexpr float resolution = resolution_micrometers / 1000.0f; // [mm]

        // TODO: Change static_cast to linear interpolation?
        const unsigned int x = static_cast<unsigned int>(x_millis / resolution) % size;
        const unsigned int y = static_cast<unsigned int>(y_millis / resolution) % size;
        const unsigned int z = static_cast<unsigned int>(z_millis / resolution) % size;

        //return matrix[x][y][z];
    }

private:
    std::array<std::array<std::array<float, size>, size>, size> matrix;
};

#endif // VOLUME_H

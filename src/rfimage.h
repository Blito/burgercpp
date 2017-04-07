#ifndef RFIMAGE_H
#define RFIMAGE_H

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <array>
#include <units/units.h>

#include "psf.h"

/**
 * Radio-frequency image.
 *
 * Stores the resulting echoes of the ultrasound colliding with tissues.
 * Each column should gather information from a single transducer element.
 * The number of rows is calculated automatically according to maximum travel
 * time of the ultrasound pulse, psf's axial resolution, and average speed of sound.
 */
template <unsigned int columns, unsigned int max_travel_time /*μs*/, unsigned int axial_resolution /*μm*/, unsigned int speed_of_sound = 1500 /*μm/μs*/>
class rf_image
{
public:
    rf_image(units::length::millimeter_t radius, units::angle::radian_t angle) :
        intensities(max_rows, columns, CV_32FC1),
        conv_axial_buffer(max_rows, columns, CV_32FC1),
        scan_converted(400, 500, CV_32FC1)
    {
        std::cout << "rf_image: " << max_rows << ", " << columns << std::endl;
        create_mapping(radius, angle, columns, max_rows);
    }

    void add_echo(const unsigned int column, const float echo, const units::time::microsecond_t micros_from_source)
    {
        const units::dimensionless::dimensionless_t row = micros_from_source / (axial_resolution_ / speed_of_sound_);

        if (row < max_rows)
        {
            intensities.at<float>(row, column) += echo;
        }
    }

    // get the delta time that represents a pixel (row resolution in time)
    constexpr units::time::microsecond_t get_dt() const
    {
        return axial_resolution / speed_of_sound_;
    }

    constexpr units::time::microsecond_t micros_traveled(units::length::micrometer_t microm_from_source) const
    {
        return microm_from_source / speed_of_sound_;
    }

    // Transforms the rf image by doing a fast approximation of the envelope function
    void envelope()
    {
        // Travel through each column looking for concave peaks.
        // Then, recalculate the points between the peaks as the linear interpolation of the absolute values of those peaks.
        // This should work as a fast approximation of the hilbert transform over the rf signal.

        for (size_t column = 0; column < columns; column++)
        {
            bool ascending = intensities.at<float>(0, column) < intensities.at<float>(1, column);
            size_t last_peak_pos = 0;
            float last_peak = intensities.at<float>(last_peak_pos, column);
            for (size_t i = 1; i < max_rows-1; i++)
            {
                if (intensities.at<float>(i, column) < intensities.at<float>(i+1, column))
                {
                    ascending = true;
                }
                else if (ascending)
                // if it was ascending and now descended, we found a concave point at i
                {
                    ascending = false;
                    const float new_peak = std::abs(intensities.at<float>(i, column));

                    // lerp last_peak -> new_peak over last_peak_pos -> i (new_peak_pos)
                    for (size_t j = last_peak_pos; j < i; j++)
                    {
                        const float alpha = (static_cast<float>(j) - static_cast<float>(last_peak_pos)) /
                                            (static_cast<float>(i) - static_cast<float>(last_peak_pos));

                        intensities.at<float>(j, column) = last_peak * (1-alpha) + new_peak * alpha;
                    }

                    last_peak_pos = i;
                    last_peak = new_peak;
                }
            }
        }
    }

    template <typename psf_>
    void convolve(const psf_ & p)
    {
        // Convolve using only axial kernel and store in intermediate buffer
        for (int col = 0; col < intensities.cols; col++) //each column is a different TE
        {
            for (int row = p.get_axial_size(); row < intensities.rows - p.get_axial_size(); row++) // each row holds information from along a ray
            {
                float convolution = 0;
                for (int kernel_i = 0; kernel_i < p.get_axial_size(); kernel_i++)
                {
                    convolution += intensities.at<float>(row + kernel_i, col) * p.axial_kernel[kernel_i];
                }
                conv_axial_buffer.at<float>(row,col) = convolution;
            }
        }

        // Convolve intermediate buffer using lateral kernel
        for (int row = p.get_axial_size(); row < conv_axial_buffer.rows - p.get_axial_size(); row++) // each row holds information from along a ray
        {
            for (int col = p.get_lateral_size() / 2; col < conv_axial_buffer.cols - p.get_lateral_size(); col++) //each column is a different TE
            {
                float convolution = 0;
                for (int kernel_i = 0; kernel_i < p.get_lateral_size(); kernel_i++)
                {
                    convolution += conv_axial_buffer.at<float>(row, col + kernel_i) * p.lateral_kernel[kernel_i];
                }
                intensities.at<float>(row,col) = convolution;
            }
        }
    }

    void postprocess()
    {
        double min, max;
        cv::minMaxLoc(intensities, &min, &max);

        save("prelog.png");

        for (size_t i = 0; i < max_rows * columns; i++)
        {
            intensities.at<float>(i) = std::log10(intensities.at<float>(i)+1)/std::log10(max+1);
        }

        // apply scan conversion using preprocessed mapping
        constexpr float invalid_color = 0.0f;
        cv::remap(intensities, scan_converted, map_y, map_x, CV_INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(invalid_color));
    }

    void save(const std::string & filename) const
    {
        cv::imwrite(filename, intensities);
    }

    void show() const
    {
        //cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE );
        //cv::imshow("Display window", image );

        cv::namedWindow("Scan Converted", cv::WINDOW_AUTOSIZE );
        cv::imshow("Scan Converted", scan_converted );

        cv::waitKey(16);
    }

    void clear()
    {
        intensities.setTo(0.0f);
    }

    void print(size_t column) const
    {
        for (size_t i = 0; i < max_rows; i++)
        {
            std::cout << intensities.at<float>(i, column) << ", ";
        }
        std::cout << std::endl;
    }

private:
    // Create constexpr variables to use type operations later
    static constexpr units::velocity::meters_per_second_t speed_of_sound_ = units::velocity::meters_per_second_t(speed_of_sound);
    static constexpr units::length::micrometer_t axial_resolution_ = units::length::micrometer_t(axial_resolution);

    static constexpr unsigned int max_rows = (speed_of_sound * max_travel_time) / axial_resolution;

    // fills map_x and map_y
    void create_mapping(units::length::millimeter_t radius, units::angle::radian_t total_angle, unsigned int rf_width, unsigned int rf_height)
    {
        // ratio to convert from mm to px
        float ratio = (max_travel_time * speed_of_sound * 0.001f + radius.to<float>() - radius.to<float>() * std::cos(total_angle.to<float>()/2.0)) / scan_converted.rows;

        // distance to transducer center going from the edge of the scan_converted image
        units::length::millimeter_t shift_y = radius * std::cos(total_angle.to<float>() / 2.0f);

        // horizontal center of the image
        float half_width = (float)scan_converted.cols / 2.0f;

        map_x.create(scan_converted.size(), CV_32FC1);
        map_y.create(scan_converted.size(), CV_32FC1);

        for( int j = 0; j < scan_converted.cols; j++ )
        {
            for( int i = 0; i < scan_converted.rows; i++ )
            {
                float fi = static_cast<float>(i)+shift_y.to<float>()/ratio;
                float fj = static_cast<float>(j)-half_width;

                // distance from the point i,j to the center of the transducer, located in -shift_y,half_width
                float r = std::sqrt(std::pow(fi,2.0f) + std::pow(fj,2.0f));

                // angle of the vector (i+shift_y, j-half_width) against the half center line
                units::angle::radian_t angle = units::angle::radian_t(std::atan2(fj, fi));

                // Invalid values are OK here. cv::remap checks for them and assigns them a special color.
                map_x.at<float>(i,j) = (r*ratio-radius.to<float>())/(max_travel_time*speed_of_sound*0.001f) * (float)rf_height;
                map_y.at<float>(i,j) = ((angle - (-total_angle/2)) / (total_angle)) * (float)rf_width;
            }
        }
    }

    cv::Mat intensities;
    cv::Mat scan_converted;
    cv::Mat conv_axial_buffer;

    // Mapping used for scan conversion
    cv::Mat map_x, map_y;
};

// http://stackoverflow.com/a/22414046
template <unsigned int columns, unsigned int max_travel_time, unsigned int axial_resolution, unsigned int speed_of_sound>
constexpr units::velocity::meters_per_second_t rf_image<columns, max_travel_time, axial_resolution, speed_of_sound>::speed_of_sound_;

template <unsigned int columns, unsigned int max_travel_time, unsigned int axial_resolution, unsigned int speed_of_sound>
constexpr units::length::micrometer_t rf_image<columns, max_travel_time, axial_resolution, speed_of_sound>::axial_resolution_;

#endif // RFIMAGE_H

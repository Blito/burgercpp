#ifndef RFIMAGE_H
#define RFIMAGE_H

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <array>

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
    rf_image() :
        image(max_rows, columns, CV_32FC1)
    {
        std::cout << "rf_image: " << max_rows << ", " << columns << std::endl;
    }

    void add_echo(const unsigned int column, const float echo_intensity, const float micros_from_source)
    {
        const auto row = static_cast<unsigned int>(static_cast<float>(micros_from_source) / (static_cast<float>(axial_resolution) / static_cast<float>(speed_of_sound)));

        if (row < max_rows)
        {
            image.at<float>(row, column) += echo_intensity;
        }
    }

    // get the delta time that represents a pixel (row resolution in time)
    constexpr float get_dt() const // [μs]
    {
        return static_cast<float>(axial_resolution)/static_cast<float>(speed_of_sound);
    }

    constexpr float micros_traveled(unsigned int microm_from_source) const
    {
        return static_cast<float>(microm_from_source)/static_cast<float>(speed_of_sound);
    }

    // Transforms the rf image by doing a fast approximation of the envelope function
    void envelope()
    {
        // Travel through each column looking for concave peaks.
        // Then, recalculate the points between the peaks as the linear interpolation of the absolute values of those peaks.
        // This should work as a fast approximation of the hilbert transform over the rf signal.

        for (size_t column = 0; column < columns; column++)
        {
            bool ascending = image.at<float>(0, column) < image.at<float>(1, column);
            size_t last_peak_pos = 0;
            float last_peak = image.at<float>(last_peak_pos, column);
            for (size_t i = 1; i < max_rows-1; i++)
            {
                if (image.at<float>(i, column) < image.at<float>(i+1, column))
                {
                    ascending = true;
                }
                else if (ascending)
                // if it was ascending and now descended, we found a concave point at i
                {
                    ascending = false;
                    const float new_peak = std::abs(image.at<float>(i, column));

                    // lerp last_peak -> new_peak over last_peak_pos -> i (new_peak_pos)
                    for (size_t j = last_peak_pos; j < i; j++)
                    {
                        const float alpha = (static_cast<float>(j) - static_cast<float>(last_peak_pos)) /
                                            (static_cast<float>(i) - static_cast<float>(last_peak_pos));

                        image.at<float>(j, column) = last_peak * (1-alpha) + new_peak * alpha;
                    }

                    last_peak_pos = i;
                    last_peak = new_peak;
                }
            }
        }
    }

    void show() const
    {
        cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE );
        cv::imshow("Display window", image );

        cv::waitKey(30);
    }

    void clear()
    {
        image.setTo(0.0f);
    }

    void print(size_t column) const
    {
        for (size_t i = 0; i < max_rows; i++)
        {
            std::cout << image.at<float>(i, column) << ", ";
        }
        std::cout << std::endl;
    }

private:
    static constexpr unsigned int max_rows = (max_travel_time * speed_of_sound) / axial_resolution;

    cv::Mat image;
};

#endif // RFIMAGE_H

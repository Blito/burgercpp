#ifndef RFIMAGE_H
#define RFIMAGE_H

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <array>

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

    // get the delta time that represents a pixel
    constexpr float get_dt() const // [μs]
    {
        return static_cast<float>(axial_resolution)/static_cast<float>(speed_of_sound);
    }

    constexpr float micros_traveled(unsigned int microm_from_source) const
    {
        return static_cast<float>(microm_from_source)/static_cast<float>(speed_of_sound);
    }


    void show()
    {
        cv::Mat normalized;
        image.convertTo(normalized, -1, 0.5, 0.5);

        cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE );
        cv::imshow("Display window", normalized );

        cv::waitKey(0);
    }

    void clear()
    {
        image.setTo(0.0f);
    }

    void print(size_t column)
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

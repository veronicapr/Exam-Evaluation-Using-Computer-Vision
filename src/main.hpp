/*!
 * @brief Detects and outputs score of a multiple choice exam or exam part from video feed.
 *
 * @authors Miguel Ferreira, Veronica Rocha
 */

#ifndef	MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <opencv2/opencv.hpp>
#include "zbar.h"

/*!
 * @brief Detector general information.
 */
namespace mc_detector
{
    /*!
     * @typedef Main window name
     */
    const char* ORIGINAL_IMAGE = "Multiple Choice Detector";
    /*!
     * @typedef Work window name
     */
    const char* WORK_IMAGE = "Work image";
    /*!
     * @typedef Work window name
     */
    const char* SCORE_IMAGE = "Score image";
    /*!
     * @typedef Work window name
     */
    const char* TABLE_DETECTED = "Table Detected";

    struct Update_Info
    {
            /*! @brief Number answer lines per column. */
            unsigned int answer_lines;
            /*! @brief Multiple choice columns. */
            unsigned int mc_columns;
            /*! @brief Number of multiple choice answers. */
            unsigned int mc_answers;
            /*! @brief True-False columns. */
            unsigned int tf_columns;
    };
}
namespace table_search
{
    const unsigned int CANNY_UPPER_THRESHOLD = 30;

}
namespace answer_table
{
    /*!
     * @struct Table main.hpp "main.hpp"
     *
     * @brief Contains information over the answer table.
     */
    struct Table_Info
    {
            /*! @brief Table image. */
            cv::Mat table;
            /*! @brief Table original corners. */
            cv::Point corners[4];
            /*! @brief Number answer lines per column. */
            unsigned int answer_lines;
            /*! @brief Multiple choice columns. */
            unsigned int mc_columns;
            /*! @brief Number of multiple choice answers. */
            unsigned int mc_answers;
            /*! @brief True-False columns. */
            unsigned int tf_columns;
            /*! @brief Table found flag. */
            bool table_found;
    };
    /*!
     * @struct Table main.hpp "main.hpp"
     *
     * @brief Contains information over the answer table.
     */
    struct Answer_Info
    {
            /*! @brief Total score. */
            double score;
            /*! @brief Multiple choice answers. */
            std::vector<std::vector<int>> multiple_choice;
            /*! @brief True false answers. */
            std::vector<std::vector<int>> true_false;
    };
}
namespace correction
{
    struct Correction_Info
    {
            /*! @brief Table image. */
            cv::Mat corrected_table;
            /*! @brief Multiple choice score. */
            double multiple_choice_value;
            /*! @brief Multiple choice score. */
            double multiple_choice_error;
            /*! @brief Multiple choice score. */
            double true_false_value;
            /*! @brief Multiple choice score. */
            double true_false_error;
            /*! @brief Multiple choice answers. */
            std::vector<int> multiple_choice;
            /*! @brief True false answers. */
            std::vector<int> true_false;
    };
}
/*!
 * @brief Detector general information.
 */
namespace mc_detector
{

    void printScore(cv::Mat score, answer_table::Answer_Info answer_info);

    void fillMainWindow(cv::Mat camera_frame);

    void fillTableWindow(cv::Mat score, cv::Mat table, cv::Mat work_image, mc_detector::Update_Info update_info);

    void trackbarAnswerLines(int pos, void *userdata);

    void trackbarMCColumns(int pos, void *userdata);

    void trackbarMCAnswers(int pos, void *userdata);

    void trackbarTFColumns(int pos, void *userdata);

    void setUpdate(mc_detector::Update_Info* update_info, answer_table::Table_Info* table_info);

    void filterImage(cv::Mat* work_image);
}
namespace table_search
{

    void findTableArea(answer_table::Table_Info* table_info, cv::Mat* original, cv::Mat* work_image);

    void getBiggestRectangle(std::vector<std::vector<cv::Point>>* rectangles, unsigned int* rectangle_index);

    void filterForRectanglePolygons(std::vector<std::vector<cv::Point>>* contours, std::vector<std::vector<cv::Point>>* rectangles);

    bool checkIfTable(answer_table::Table_Info* table_info, cv::Mat* work_image, std::vector<cv::Point>* rectangle);

    void deformTable(cv::Mat* work_image, cv::Mat* table_cut, std::vector<cv::Point>* rectangle);
}
namespace answer_table
{

    void baseTableDimensions(answer_table::Table_Info* table_info);

    void redoAnswerTable(answer_table::Table_Info* table_info, answer_table::Answer_Info* answer_info);

    void getChoices(answer_table::Answer_Info* answer_info, answer_table::Table_Info* table_info);

    bool getAnswers(cv::Mat* table, cv::Rect cut_rectangle);
}
namespace correction
{

    void baseCorrectionInfo(correction::Correction_Info* correction_info);

    void correctAnswers(correction::Correction_Info* correction_info, answer_table::Answer_Info* answer_info, answer_table::Table_Info* table_info);

    void wrongAnswer(cv::Mat* table, cv::Rect cut_rectangle);

    void correctAnswer(cv::Mat* table, cv::Rect cut_rectangle);
}

#endif

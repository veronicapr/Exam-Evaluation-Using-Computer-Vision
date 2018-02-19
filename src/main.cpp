/*!
 * @brief Detects and outputs score of a multiple choice exam or exam part from video feed.
 *
 * @authors Miguel Ferreira, Veronica Rocha
 */

#include "main.hpp"

unsigned int slider_answer_lines;
unsigned int slider_mc_columns;
unsigned int slider_mc_answers;
unsigned int slider_tf_columns;

/* ========================= ========================= Main ========================= ========================= */
/*!
 * @brief Program main
 */
int main(int argc, char* argv[]) {
    // computer camera feed
    cv::VideoCapture camera_feed;
    // video frame to be processed and image used for the process
    cv::Mat camera_frame, work_image, table_detected, score;
    // update info
    mc_detector::Update_Info update_info = { };
    // table information structure
    answer_table::Table_Info table_info = { };
    // answer information structure
    answer_table::Answer_Info answer_info = { };
    // correction information structure
    correction::Correction_Info correction_info = { };
    // opens a camera feed. Stops program is it fails to open one
    if (!camera_feed.open(0)) {
        std::cerr << "Unable to open video capture" << std::endl;
        return -1;
    }
    // initialises table information
    answer_table::baseTableDimensions(&table_info);
    // initialises correction information
    correction::baseCorrectionInfo(&correction_info);
    // set update info according to table info
    mc_detector::setUpdate(&update_info, &table_info);
    // processes camera feed until "escape" is pressed.
    while (true) {
        // reset table to table dimensions
        answer_table::redoAnswerTable(&table_info, &answer_info);
        // reset score
        score = cv::Mat::zeros(150, 500, CV_8UC3);
        // reset score
        table_detected = cv::Mat::zeros(100, 200, CV_8UC3);
        // create windows
        cv::namedWindow(mc_detector::ORIGINAL_IMAGE, CV_WINDOW_AUTOSIZE);
        cv::namedWindow(mc_detector::WORK_IMAGE, CV_WINDOW_AUTOSIZE);
        cv::namedWindow(mc_detector::SCORE_IMAGE, CV_WINDOW_AUTOSIZE);
        cv::namedWindow(mc_detector::TABLE_DETECTED, CV_WINDOW_AUTOSIZE);
        // gets new frame
        camera_feed >> camera_frame;
        // frame is empty, move on
        if (camera_frame.empty()) continue;
        // convert image to grey scale
        cv::cvtColor(camera_frame, work_image, cv::COLOR_BGR2GRAY);
        // apply filter to image
        mc_detector::filterImage(&work_image);
        // find table area
        table_search::findTableArea(&table_info, &camera_frame, &work_image);
        // if table is found
        if (table_info.table_found) {
            // trade work image for table
            correction_info.corrected_table = table_info.table;
            // change to colour BGR
            cv::cvtColor(correction_info.corrected_table, correction_info.corrected_table, cv::COLOR_GRAY2BGR);
            // get answers
            answer_table::getChoices(&answer_info, &table_info);
            // correct answers
            correction::correctAnswers(&correction_info, &answer_info, &table_info);
            // table corrected
            table_detected = correction_info.corrected_table;
        }
        // update table info
        table_info.answer_lines = slider_answer_lines;
        table_info.mc_answers = slider_mc_answers;
        table_info.mc_columns = slider_mc_columns;
        table_info.tf_columns = slider_tf_columns;
        // print score
        mc_detector::printScore(score, answer_info);
        // fill main window
        mc_detector::fillMainWindow(camera_frame);
        // fill table window
        mc_detector::fillTableWindow(score, table_detected, work_image, update_info);
        // awaits program termination
        if (cv::waitKey(25) == 27) break;
    }
}

/* ========================= ========================= MCDetector ========================= ========================= */
void mc_detector::printScore(cv::Mat score, answer_table::Answer_Info answer_info) {
    char grade[10];
    sprintf(grade, "%2.3f", answer_info.score);
    cv::putText(score, grade, cv::Point(0, 100), cv::FONT_HERSHEY_COMPLEX, 3, cv::Scalar(255, 255, 255, 255), 1);
}

void mc_detector::fillMainWindow(cv::Mat camera_frame) {
    cv::imshow(mc_detector::ORIGINAL_IMAGE, camera_frame);
}

void mc_detector::fillTableWindow(cv::Mat score, cv::Mat table, cv::Mat work_image, mc_detector::Update_Info update_info) {
    cv::createTrackbar("Answer Lines", mc_detector::WORK_IMAGE, (int*) &slider_answer_lines, 15, mc_detector::trackbarAnswerLines, &update_info);
    cv::setTrackbarMin("Answer Lines", mc_detector::WORK_IMAGE, 3);
    cv::createTrackbar("MC Columns", mc_detector::WORK_IMAGE, (int*) &slider_mc_columns, 3, mc_detector::trackbarMCColumns, &update_info);
    cv::createTrackbar("MC Answers", mc_detector::WORK_IMAGE, (int*) &slider_mc_answers, 5, mc_detector::trackbarMCAnswers, &update_info);
    cv::createTrackbar("TF Columns", mc_detector::WORK_IMAGE, (int*) &slider_tf_columns, 3, mc_detector::trackbarTFColumns, &update_info);
    cv::imshow(mc_detector::WORK_IMAGE, work_image);
    cv::imshow(mc_detector::SCORE_IMAGE, score);
    cv::imshow(mc_detector::TABLE_DETECTED, table);
}

void mc_detector::trackbarAnswerLines(int pos, void *userdata) {
    mc_detector::Update_Info* update_info = (mc_detector::Update_Info*) userdata;
    update_info->answer_lines = slider_answer_lines;
}

void mc_detector::trackbarMCColumns(int pos, void *userdata) {
    mc_detector::Update_Info* update_info = (mc_detector::Update_Info*) userdata;
    update_info->mc_columns = slider_mc_columns;
}

void mc_detector::trackbarMCAnswers(int pos, void *userdata) {
    mc_detector::Update_Info* update_info = (mc_detector::Update_Info*) userdata;
    update_info->mc_answers = slider_mc_answers;
}

void mc_detector::trackbarTFColumns(int pos, void *userdata) {
    mc_detector::Update_Info* update_info = (mc_detector::Update_Info*) userdata;
    update_info->tf_columns = slider_tf_columns;
}

void mc_detector::setUpdate(mc_detector::Update_Info* update_info, answer_table::Table_Info* table_info) {
    update_info->answer_lines = table_info->answer_lines;
    update_info->mc_answers = table_info->mc_answers;
    update_info->mc_columns = table_info->mc_columns;
    update_info->tf_columns = table_info->tf_columns;
    slider_answer_lines = table_info->answer_lines;
    slider_mc_answers = table_info->mc_answers;
    slider_mc_columns = table_info->mc_columns;
    slider_tf_columns = table_info->tf_columns;
}

void mc_detector::filterImage(cv::Mat* work_image) {
    // blurred image
    cv::Mat blurred_image;
    // threshold result
    cv::Mat threshold;
    // apply blur
    cv::GaussianBlur(*work_image, blurred_image, cv::Size(3, 3), 0.0);
    // apply threshold
    cv::threshold(blurred_image, threshold, 0.0, 255.0, cv::THRESH_OTSU);
    // copy threshold to grey_image
    threshold.copyTo(*work_image);
}
/* ========================= ========================= Table Search ========================= ========================= */
void table_search::findTableArea(answer_table::Table_Info* table_info, cv::Mat* original, cv::Mat* work_image) {
    // image with edge detection
    cv::Mat edge_detection;
    // found contours
    std::vector<std::vector<cv::Point>> contours;
    // found rectangles
    std::vector<std::vector<cv::Point>> rectangles;
    // approximate polygon from contours/current rectangle
    std::vector<cv::Point> rectangle;
    // max rectangle index
    unsigned int biggest_index;
    // force table found flag to false
    table_info->table_found = false;
    // set edge detection image to same size as work_image
    edge_detection = cv::Mat(work_image->size(), CV_8U);
    // apply edge detection to the image
    cv::Canny(*work_image, edge_detection, 0, table_search::CANNY_UPPER_THRESHOLD, 5);
    // dilate contours
    cv::dilate(edge_detection, edge_detection, cv::Mat());
    // find contours on the image
    cv::findContours(edge_detection, contours, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);
    // filter for rectangles
    table_search::filterForRectanglePolygons(&contours, &rectangles);
    // draw rectangles
    cv::drawContours(*original, rectangles, -1, cv::Scalar(0, 255, 0, 255), 1);
    // check if any rectangles were found, exits if none
    if (rectangles.empty()) return;
    // cycle trough all rectangles until table is found
    while (!rectangles.empty()) {
        // get biggest rectangle index
        table_search::getBiggestRectangle(&rectangles, &biggest_index);
        // get first polygon
        rectangle = rectangles[biggest_index];
        // check if biggest rectangle is the table
        if (table_search::checkIfTable(table_info, work_image, &rectangle)) {
            // store the found rectangle
            table_info->corners[0] = rectangle[0];
            table_info->corners[1] = rectangle[1];
            table_info->corners[2] = rectangle[2];
            table_info->corners[3] = rectangle[3];
            // set flag to true
            table_info->table_found = true;
            return;
        }
        // delete biggest entry if it isn't the table
        rectangles.erase(rectangles.begin() + biggest_index);
    }
}

void table_search::filterForRectanglePolygons(std::vector<std::vector<cv::Point>>* contours, std::vector<std::vector<cv::Point>>* rectangles) {
    // approximate polygon from contours
    std::vector<cv::Point> polygon;
    // points for cosine calculations
    cv::Point a, b, c;
    // maximum cosine, and current cosine of the angle between joint edges
    float max_cosine, cosine;
    // length variables for the sides in angle
    float dx1, dx2, dy1, dy2;
    // cycle to each found contours
    for (unsigned int i = 0; i < contours->size(); i++) {
        // determine approximate polygon according to contours
        cv::approxPolyDP(cv::Mat((*contours)[i]), polygon, cv::arcLength(cv::Mat((*contours)[i]), true) * 0.02, true);
        // filter for rectangles only
        if (polygon.size() == 4 && fabs(cv::contourArea(cv::Mat(polygon))) > 1000 && cv::isContourConvex(cv::Mat(polygon))) {
            // set max cosine as 0
            max_cosine = 0;
            // find the maximum cosine of the angle between joint edges
            for (unsigned int j = 2; j < 5; j++) {
                // get points
                a = polygon[j - 1];
                b = polygon[j % 4];
                c = polygon[j - 2];
                // calculate distance
                dx1 = b.x - a.x;
                dy1 = b.y - a.y;
                dx2 = c.x - a.x;
                dy2 = c.y - a.y;
                // calculate cosine
                cosine = (dx1 * dx2 + dy1 * dy2) / sqrt((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2) + 1e-10);
                // set max cosine
                max_cosine = MAX(max_cosine, cosine);
            }
            if (max_cosine <= 0.5) {
                rectangles->push_back(polygon);
            }
        }
    }
}

void table_search::getBiggestRectangle(std::vector<std::vector<cv::Point>>* rectangles, unsigned int* rectangle_index) {
    // approximate polygon from contours/current rectangle
    std::vector<cv::Point> polygon;
    // rectangle areas
    double rectangle_area, max_rectangle_area;
    // max rectangle index
    unsigned int max_area_index;
    // reset max area index
    max_area_index = 0;
    // get first polygon
    polygon = (*rectangles)[max_area_index];
    // calculate first polygon area
    max_rectangle_area = cv::contourArea(polygon);
    // cycle trough all others to get max area
    for (unsigned int i = 1; i < rectangles->size(); i++) {
        // get polygon i
        polygon = (*rectangles)[i];
        // calculate area
        rectangle_area = cv::contourArea(polygon);
        // check if new is bigger
        if (rectangle_area > max_rectangle_area) {
            // update
            max_area_index = i;
            max_rectangle_area = rectangle_area;
        }
    }
    // store max area index
    *rectangle_index = max_area_index;
}

bool table_search::checkIfTable(answer_table::Table_Info* table_info, cv::Mat* work_image, std::vector<cv::Point>* rectangle) {
    // table cut
    cv::Mat table_cut;
    // top corner cut
    cv::Mat top_corner;
    // top corner rectangle
    cv::Rect top_corner_rectangle;
    // temporary point (used for rotation purposes)
    cv::Point tmp;
    // pointer value
    std::vector<cv::Point> quad;
    // table rows, columns and middle values for zero percentage
    unsigned int total_rows, total_columns, total_pixels, zero_count;
    // percentage of zeros
    double zero_percentage;
    // assign rectangle
    quad = *rectangle;
    // calculate total rows and total columns of answer table
    total_rows = table_info->answer_lines + 1;
    total_columns = (table_info->mc_answers + 1) * table_info->mc_columns + 3 * table_info->tf_columns;
    // for each orientation
    for (unsigned int i = 1; i < 4; i++) {
        // rotation needed
        if (i != 0) {
            // rotate 90 degrees
            tmp = quad[0];
            quad[0] = quad[1];
            quad[1] = quad[2];
            quad[2] = quad[3];
            quad[3] = tmp;
        }
        // get a straight table
        table_search::deformTable(work_image, &table_cut, &quad);
        // determine rectangle for top corner form table dimensions
        top_corner_rectangle = cv::Rect(0, 0, table_cut.cols / total_columns, table_cut.rows / total_rows);
        // get top corner
        top_corner = table_cut(top_corner_rectangle);
        // total pixels in top_corner
        total_pixels = top_corner_rectangle.size().width * top_corner_rectangle.size().height;
        // count zeros
        zero_count = cv::countNonZero(top_corner);
        // calculate zero percentage
        zero_percentage = 1.0 - ((double) zero_count / (double) total_pixels);
        // its correct table and orientation if total black pixels over 90 percent of the corner
        if (zero_percentage > 0.8) {
            // assign current rotation to the
            *rectangle = quad;
            // store current image cut
            table_cut.copyTo(table_info->table);
            // return successful
            return true;
        }
    }
    // nothing found
    return false;
}

void table_search::deformTable(cv::Mat* work_image, cv::Mat* table_cut, std::vector<cv::Point>* rectangle) {
    // rectangle cut
    cv::Mat cut;
    // image transformation matrix
    cv::Mat transformation_matrix;
    // pointer value
    std::vector<cv::Point2f> quad;
    // wanted orientation
    std::vector<cv::Point2f> straight_quad;
    // rectangle weight and height
    double distance, weight, height;
    // assign rectangle
    quad.push_back(cv::Point2f((*rectangle)[0]));
    quad.push_back(cv::Point2f((*rectangle)[1]));
    quad.push_back(cv::Point2f((*rectangle)[2]));
    quad.push_back(cv::Point2f((*rectangle)[3]));
    // get the weight form the distance between point 0-1
    weight = cv::norm((cv::Point2f) quad[1] - (cv::Point2f) quad[0]);
    // get the weight form the distance between point 1-2
    height = cv::norm((cv::Point2f) quad[2] - (cv::Point2f) quad[1]);
    // weight as biggest distance between points 0-1 and 2-3
    distance = cv::norm((cv::Point2f) quad[3] - (cv::Point2f) quad[2]);
    if (weight < distance) weight = distance;
    // weight as biggest distance between points 1-2 and 3-0
    distance = cv::norm((cv::Point2f) quad[0] - (cv::Point2f) quad[3]);
    if (height < distance) height = distance;
    // Set table cut image size to table size
    cut = cv::Mat::zeros(weight, height, CV_8UC3);
    // straight points
    straight_quad.push_back(cv::Point2f(0, 0));
    straight_quad.push_back(cv::Point2f(0, weight));
    straight_quad.push_back(cv::Point2f(height, weight));
    straight_quad.push_back(cv::Point2f(height, 0));
    // create transformation matrix
    transformation_matrix = cv::getPerspectiveTransform(quad, straight_quad);
    // apply transformation matrix
    cv::warpPerspective(*work_image, cut, transformation_matrix, cut.size(), cv::INTER_NEAREST);
    // copy cut to table cut
    cut.copyTo(*table_cut);
}
/* ========================= ========================= Answer Table ========================= ========================= */

void answer_table::baseTableDimensions(answer_table::Table_Info* table_info) {
    // initialise corners
    table_info->corners[0] = cv::Point(0.0, 0.0);
    table_info->corners[1] = cv::Point(0.0, 0.0);
    table_info->corners[2] = cv::Point(0.0, 0.0);
    table_info->corners[3] = cv::Point(0.0, 0.0);
    // base dimensions
    table_info->answer_lines = 7;
    table_info->mc_columns = 2;
    table_info->mc_answers = 4;
    table_info->tf_columns = 1;
}

void answer_table::redoAnswerTable(answer_table::Table_Info* table_info, answer_table::Answer_Info* answer_info) {
    // temporary assignment
    std::vector<int> empty_vector;
    // size
    unsigned int size;
    // reset score
    answer_info->score = 0.0;
    // delete previous answers
    answer_info->multiple_choice.clear();
    answer_info->true_false.clear();
    // allocate space for all the questions
    size = table_info->answer_lines * table_info->mc_columns;
    for (unsigned int i = 0; i < size; i++) {
        answer_info->multiple_choice.push_back(empty_vector);
    }
    size = table_info->answer_lines * table_info->tf_columns;
    for (unsigned int i = 0; i < size; i++) {
        answer_info->true_false.push_back(empty_vector);
    }
}

void answer_table::getChoices(answer_table::Answer_Info* answer_info, answer_table::Table_Info* table_info) {
    // choice vector
    std::vector<int> current_answer;
    // section rows, columns, cut rows, columns
    unsigned int total_rows, total_columns, cut_rows, cut_columns;
    // calculate cut rows and columns
    total_rows = table_info->answer_lines + 1;
    total_columns = (table_info->mc_answers + 1) * table_info->mc_columns + 3 * table_info->tf_columns;
    cut_rows = table_info->table.rows / total_rows;
    cut_columns = table_info->table.cols / total_columns;
    // cycle trough all answers
    for (unsigned int col = 0; col < total_columns; col++) {
        // verify if on multiple choices or true and false
        if (col < ((table_info->mc_answers + 1) * table_info->mc_columns)) {
            // jump multiples of number o answers + 1
            if ((col % (table_info->mc_answers + 1)) == 0) {
                continue;
            }
            // get multiple choice answer
            for (unsigned int row = 1; row < total_rows; row++) {
                // check answer
                if (answer_table::getAnswers(&(table_info->table), cv::Rect(col * cut_columns, row * cut_rows, cut_columns, cut_rows))) {
                    // store answer
                    current_answer = answer_info->multiple_choice[row - 1 + (table_info->answer_lines * (col / (table_info->mc_answers + 1)))];
                    current_answer.push_back(col % (table_info->mc_answers + 1));
                    answer_info->multiple_choice[row - 1 + (table_info->answer_lines * (col / (table_info->mc_answers + 1)))] = current_answer;
                }
            }
        } else {
            // jump multiples of 3
            if (((col - ((table_info->mc_answers + 1) * table_info->mc_columns)) % 3) == 0) {
                continue;
            }
            // get multiple choice answer
            for (unsigned int row = 1; row < total_rows; row++) {
                // check answer
                if (answer_table::getAnswers(&(table_info->table), cv::Rect(col * cut_columns, row * cut_rows, cut_columns, cut_rows))) {
                    // store answer
                    current_answer = answer_info->true_false[row - 1 + (table_info->answer_lines * ((col - ((table_info->mc_answers + 1) * table_info->mc_columns)) / 3))];
                    current_answer.push_back((col - ((table_info->mc_answers + 1) * table_info->mc_columns)) % 3);
                    answer_info->true_false[row - 1 + (table_info->answer_lines * ((col - ((table_info->mc_answers + 1) * table_info->mc_columns)) / 3))] = current_answer;
                }
            }
        }
    }
}

bool answer_table::getAnswers(cv::Mat* table, cv::Rect cut_rectangle) {
    // rectangle
    cv::Rect smaller_area;
    // top corner cut
    cv::Mat cut;
    // pixel count
    unsigned int total_pixels, zero_count;
    // percentage of zeros
    double zero_percentage;
    // reduce area
    smaller_area = cv::Rect(cut_rectangle.x, cut_rectangle.y, cut_rectangle.size().width * 0.85, cut_rectangle.size().height * 0.85);
    // get top corner
    cut = (*table)(smaller_area);
    // total pixels in top_corner
    total_pixels = smaller_area.size().width * smaller_area.size().height;
    // count zeros
    zero_count = cv::countNonZero(cut);
    // calculate zero percentage
    zero_percentage = 1.0 - ((double) zero_count / (double) total_pixels);
    // check if there is an answer
    if ((zero_percentage > 0.09) && (zero_percentage < 0.20)) return true;
    else return false;
}
/* ========================= ========================= Correction ========================= ========================= */

void correction::baseCorrectionInfo(correction::Correction_Info* correction_info) {
    correction_info->multiple_choice_value = 20 * 0.6 / 14;
    correction_info->multiple_choice_error = (20 * 0.6 / 14) / 3;
    correction_info->true_false_value = 20 * 0.4 / 7;
    correction_info->true_false_error = (20 * 0.4 / 7) / 3;
    correction_info->multiple_choice = {3, 4, 1, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1};
    correction_info->true_false = {1, 1, 1, 1, 1, 1, 1};
}

void correction::correctAnswers(correction::Correction_Info* correction_info, answer_table::Answer_Info* answer_info, answer_table::Table_Info* table_info) {
    // choice vector
    std::vector<int> current_answer;
    // total questions for cycle, size of answer array
    unsigned int total_questions, size;
    // cut rows, columns, rectangle corners
    unsigned int cut_rows, cut_columns, rect_x, rect_y;
    // calculate cut rows and columns
    cut_rows = correction_info->corrected_table.rows / (table_info->answer_lines + 1);
    cut_columns = correction_info->corrected_table.cols / ((table_info->mc_answers + 1) * table_info->mc_columns + 3 * table_info->tf_columns);
    // set total questions for multiple choices
    total_questions = table_info->answer_lines * table_info->mc_columns;
    for (unsigned int question = 0; question < total_questions; question++) {
        // get question vector
        current_answer = answer_info->multiple_choice[question];
        // check answer
        if (current_answer.size() == 1) {
            rect_x = (current_answer[0] * cut_columns) + ((question / 7) * (table_info->mc_answers + 1) * cut_columns);
            rect_y = ((1 + (question % 7)) * cut_rows);
            if (current_answer[0] == correction_info->multiple_choice[question]) {
                answer_info->score += correction_info->multiple_choice_value;
                correction::correctAnswer(&(correction_info->corrected_table), cv::Rect(rect_x, rect_y, cut_columns, cut_rows));
            } else {
                answer_info->score -= correction_info->multiple_choice_error;
                correction::wrongAnswer(&(correction_info->corrected_table), cv::Rect(rect_x, rect_y, cut_columns, cut_rows));
            }
        } else if (current_answer.size() != 0) {
            answer_info->score -= correction_info->multiple_choice_error;
            for (unsigned int i = 0; i < current_answer.size(); i++) {
                rect_x = (current_answer[i] * cut_columns) + ((question / 7) * (table_info->mc_answers + 1) * cut_columns);
                rect_y = ((1 + (question % 7)) * cut_rows);
                correction::wrongAnswer(&(correction_info->corrected_table), cv::Rect(rect_x, rect_y, cut_columns, cut_rows));
            }
        }
    }
    // set total questions for multiple choices
    total_questions = table_info->answer_lines * table_info->tf_columns;
    for (unsigned int question = 0; question < total_questions; question++) {
        // get question vector
        current_answer = answer_info->true_false[question];
        // check answer
        if (current_answer.size() == 1) {
            rect_x = (current_answer[0] * cut_columns) + ((question / 7) * 3 * cut_columns) + (cut_columns * table_info->mc_columns * (table_info->mc_answers + 1));
            rect_y = ((1 + (question % 7)) * cut_rows);
            if (current_answer[0] == correction_info->true_false[question]) {
                answer_info->score += correction_info->true_false_value;
                correction::correctAnswer(&(correction_info->corrected_table), cv::Rect(rect_x, rect_y, cut_columns, cut_rows));
            } else {
                answer_info->score -= correction_info->true_false_error;
                correction::wrongAnswer(&(correction_info->corrected_table), cv::Rect(rect_x, rect_y, cut_columns, cut_rows));
            }
        } else if (current_answer.size() != 0) {
            answer_info->score -= correction_info->true_false_error;
            for (unsigned int i = 0; i < current_answer.size(); i++) {
                rect_x = (current_answer[i] * cut_columns) + ((question / 7) * 3 * cut_columns) + (cut_columns * table_info->mc_columns * (table_info->mc_answers + 1));
                rect_y = ((1 + (question % 7)) * cut_rows);
                correction::wrongAnswer(&(correction_info->corrected_table), cv::Rect(rect_x, rect_y, cut_columns, cut_rows));
            }
        }
    }
}

void correction::wrongAnswer(cv::Mat* table, cv::Rect cut_rectangle) {
    // draw first line
    cv::line(*table, cv::Point(cut_rectangle.x + cut_rectangle.width / 4, cut_rectangle.y + cut_rectangle.height / 4),
        cv::Point(cut_rectangle.x + 3 * cut_rectangle.width / 4, cut_rectangle.y + 3 * cut_rectangle.height / 4), cv::Scalar(0, 0, 255, 255), 2);
    // draw second line
    cv::line(*table, cv::Point(cut_rectangle.x + cut_rectangle.width / 4, cut_rectangle.y + cut_rectangle.height / 4),
        cv::Point(cut_rectangle.x + 3 * cut_rectangle.width / 4, cut_rectangle.y + cut_rectangle.height / 4), cv::Scalar(0, 0, 255, 255), 2);
}

void correction::correctAnswer(cv::Mat* table, cv::Rect cut_rectangle) {
    // draw circle
    cv::circle(*table, cv::Point(cut_rectangle.x + cut_rectangle.width / 2, cut_rectangle.y + cut_rectangle.height / 2), cut_rectangle.height / 4, cv::Scalar(0, 255, 0, 255), 2);
}

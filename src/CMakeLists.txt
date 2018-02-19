cmake_minimum_required(VERSION 3.9)

project( Multiple_Choice_Detector )

find_package( OpenCV REQUIRED )

add_executable( MCDetector main.cpp main.hpp )

target_link_libraries( MCDetector ${OpenCV_LIBS} ${ZBAR_LIBRARIES} zbar )
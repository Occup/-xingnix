cmake_minimum_required( VERSION 2.8 )
project( Xingnix )

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )
# 添加cmake模块以使用ceres库
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

find_package( Ceres REQUIRED )
include_directories( ${CERES_INCLUDE_DIRS} )

find_package( G2O REQUIRED )
include_directories( ${G2O_INCLUDE_DIRS} "/usr/include/eigen3")

find_package( OpenCV 3.1 REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(EigenQuaternion EigenQuaternion.cpp)

add_executable( AccelAnalyse AccelAnalyse.cpp )
target_link_libraries( AccelAnalyse ${CERES_LIBRARIES} )

add_executable( GearAnalyse GearAnalyse.cpp )
target_link_libraries( GearAnalyse g2o_core g2o_stuff)

add_executable( ImageManipulate ImageManipulate.cpp )
target_link_libraries( ImageManipulate ${OpenCV_LIBS} )

add_executable( FeatureExtraction FeatureExtraction.cpp  )
target_link_libraries( FeatureExtraction ${OpenCV_LIBS} )

add_executable( PoseEstimation PoseEstimation.cpp )
target_link_libraries( PoseEstimation ${OpenCV_LIBS} )

add_executable( Triangulate Triangulate.cpp )
target_link_libraries( Triangulate ${OpenCV_LIBS} )

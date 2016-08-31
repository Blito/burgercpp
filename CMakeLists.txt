project(burgercpp)

cmake_minimum_required(VERSION 3.2.0)

# C++14 support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# Sources
set (burgercpp_SRC
        "src/main.cpp")

add_executable(burgercpp ${burgercpp_SRC})

# Bullet includes
set(BULLET_SRC_DIR "" CACHE PATH "")
set(BULLET_LIB_DIR "" CACHE PATH "")
message("Bullet source directory: ${BULLET_SRC_DIR}")
message("Bullet lib directory: ${BULLET_LIB_DIR}")
include_directories(${BULLET_SRC_DIR})

target_link_libraries(burgercpp
    ${BULLET_LIB_DIR}/libBulletDynamics.a
    ${BULLET_LIB_DIR}/libBulletCollision.a
    ${BULLET_LIB_DIR}/libLinearMath.a
    ${BULLET_LIB_DIR}/libBulletSoftBody.a
    )

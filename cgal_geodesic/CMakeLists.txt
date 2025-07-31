cmake_minimum_required(VERSION 3.10)
project(ExactGeodesic)

find_package(CGAL REQUIRED)

add_executable(geodesic_example geodesic_example.cpp)
target_link_libraries(geodesic_example CGAL::CGAL)
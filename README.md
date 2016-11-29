# Burgercpp

Ultrasound simulation using raytracing over surfaces.

Implementation of [Real-Time GPU-Based Ultrasound Simulation Using Deformable Mesh Models](http://sci-hub.cc/10.1109/tmi.2012.2234474) in C++ (in CPU).

---

## Example scene

![Sphere scene 3D](http://i.imgur.com/b8Acetb.png)
Test

![Sphere scene](http://i.imgur.com/OhOPUp3.png)

## Prerequisites
- C++14 compiler (tested using MinGW-w64)
- CMake 3.2

## To run
    burgercpp examples/sphere/sphere.scene

## Third Party Libraries
- [Bulletphysics](https://github.com/bulletphysics/bullet3)
- [nholthaus/units](https://github.com/nholthaus/units)
- [nlohmann/json](https://github.com/nlohmann/json)

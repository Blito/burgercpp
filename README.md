# Burgercpp

Ultrasound simulation using raytracing over surfaces.

Implementation of [Real-Time GPU-Based Ultrasound Simulation Using Deformable Mesh Models](http://sci-hub.cc/10.1109/tmi.2012.2234474) in C++ (in CPU).

---

## Some example scenes

![Sphere scene 3D](http://i.imgur.com/U71zcBx.png)
Test scene using a sphere of bone tissue inside a cube of liver tissue.

![3D Phantom reconstruction](http://i.imgur.com/zmAoiSv.png)
Simulation of a vessel structure between two layers of different tissue.

![Real patient simulation](http://i.imgur.com/QBoGAID.png)
Simulation of rays going through [Morison's pouch](https://en.wikipedia.org/wiki/Hepatorenal_recess_of_subhepatic_space). Visualization done using [this repo](https://github.com/Blito/bullet3).

---

## Prerequisites
- C++14 compiler (tested using MinGW-w64)
- CMake 3.2

## To run
    burgercpp examples/sphere/sphere.scene

## Third Party Libraries
- [Bulletphysics](https://github.com/bulletphysics/bullet3) (not included)
- [nholthaus/units](https://github.com/nholthaus/units) (header-only, included)
- [nlohmann/json](https://github.com/nlohmann/json) (header-only, included)

#ifndef SCENE_H
#define SCENE_H

#include "btBulletDynamicsCommon.h"

#include "ray.h"
#include "mesh.h"
#include "transducer.h"

#include <ctime>
#include <memory>
#include <unordered_map>
#include <array>
#include <vector>

#include <nlohmann/json.hpp>
#include <units/units.h>

class scene
{
    using transducer_ = transducer<512>;

public:
    explicit scene(const nlohmann::json & config, transducer_ & transducer);
    ~scene();

    void init();

    template<unsigned int ray_count>
    std::array<std::vector<ray_physics::segment>, ray_count> cast_rays();

    void step(float delta_time);

    units::length::millimeter_t distance(const btVector3 & from, const btVector3 & to) const;

protected:
    std::string working_dir;

    std::unordered_map<std::string, material> materials;
    std::string starting_material;

    std::vector<mesh> meshes;

    transducer_ & transducer;

    const float intensity_epsilon { 1e-8 };
    const float initial_intensity { 1.0f };

    std::array<float,3> spacing;
    std::array<float,3> origin;
    float scaling;

    void parse_config(const nlohmann::json & config);

    void create_empty_world();
    void destroy_world();

    units::length::millimeter_t distance_in_mm(const btVector3 & v1, const btVector3 & v2) const;
    btVector3 enlarge(const btVector3 & versor, float mm) const;

    class btRigidBody * add_rigidbody_from_obj(const std::string & fileName, std::array<float, 3> deltas, float scaling);

    btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;
    std::unique_ptr<btBroadphaseInterface> m_broadphase;
    std::unique_ptr<btCollisionDispatcher> m_dispatcher;
    std::unique_ptr<btConstraintSolver> m_solver;
    std::unique_ptr<btDefaultCollisionConfiguration> m_collisionConfiguration;
    std::unique_ptr<btDiscreteDynamicsWorld> m_dynamicsWorld;

    btVector3 transducer_pos;
    std::array<units::angle::degree_t, 3> transducer_dir;

    clock_t frame_start;
};

#endif // SCENE_H

#ifndef SCENE_H
#define SCENE_H

#include "btBulletDynamicsCommon.h"

#include "material.h"
#include "ray.h"

#include <ctime>
#include <memory>
#include <unordered_map>
#include <array>
#include <vector>

class scene
{
public:
    scene();
    ~scene();

    void init();

    void set_transducer(const btVector3 & position, const btVector3 & direction);

    template<unsigned int ray_count>
    std::array<std::vector<ray_physics::segment>, ray_count> cast_rays();

    void step(float delta_time);

    float distance(const btVector3 & from, const btVector3 & to);

protected:
    struct organ_properties
    {
        organ_properties(material_id mat_id)
            : mat_id(mat_id)
        {}

        material_id mat_id;
    };

    const std::string working_dir { "data/ultrasound/" };

    std::unordered_map<material_id, material, EnumClassHash> materials
    {
        { material_id::GEL, {1.99f, 1e-8, 0.0f, 0.0f, 0.0f} },
        { material_id::AIR, {0.0004f, 1.64f, 0.78f, 0.56f, 0.1f} },
        { material_id::FAT, {1.38f, 0.63f, 0.5f, 0.5f, 0.0f} },
        { material_id::BONE, {7.8f, 5.0f, 0.78f, 0.56f, 0.1f} },
        { material_id::BLOOD, {1.61f, 0.18f, 0.001f, 0.0f, 0.01f} },
        { material_id::VESSEL, {1.99f, 1.09f, 0.2f, 0.1f, 0.2f} },
        { material_id::LIVER, {1.65f, 0.7f, 0.19f, 1.0f, 0.24f} },
        { material_id::KIDNEY, {1.62f, 1.0f, 0.4f, 0.6f, 0.3f} },
        { material_id::SUPRARRENAL, {1.62f, 1.0f, 0.4f, 0.6f, 0.3f} }, // todo
        { material_id::GALLBLADDER, {1.62f, 1.0f, 0.4f, 0.6f, 0.3f} }, // todo
        { material_id::SKIN, {1.99f, 0.6f, 0.5f, 0.2f, 0.5f} },
    };

    const float frequency { 5.0f };
    const float intensity_epsilon { 1e-8 };
    const float initial_intensity { 1.0f };

    const std::array<float,3> spacing {{ 1.0f, 1.0f, 1.0f }};

    void create_empty_world();
    void destroy_world();

    float distance_in_mm(const btVector3 & v1, const btVector3 & v2) const;
    btVector3 enlarge(const btVector3 & versor, float mm) const;

    class btRigidBody * add_rigidbody_from_obj(const std::string & fileName, std::array<float, 3> deltas, float scaling);

    btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;
    std::unique_ptr<btBroadphaseInterface> m_broadphase;
    std::unique_ptr<btCollisionDispatcher> m_dispatcher;
    std::unique_ptr<btConstraintSolver> m_solver;
    std::unique_ptr<btDefaultCollisionConfiguration> m_collisionConfiguration;
    std::unique_ptr<btDiscreteDynamicsWorld> m_dynamicsWorld;

    btVector3 transducer_pos, transducer_dir;
    clock_t frame_start;
};

#endif // SCENE_H

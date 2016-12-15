#include "scene.h"

#include "objloader.h"
#include "ray.h"

#include <units/units.h>

#include <cmath>
#include <iostream>
#include <cassert>
#include <exception>

template std::array<std::vector<ray_physics::segment>, 256> scene::cast_rays<256>();

scene::scene(const nlohmann::json & config)
{
    try
    {
        parse_config(config);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error{ "Error while loading scene: " + std::string{ex.what()} };
    }

    create_empty_world();
}

scene::~scene()
{
    destroy_world();
}

void scene::init()
{
    for (auto & mesh : meshes)
    {
        const auto full_path = working_dir + mesh.filename;

        auto object = add_rigidbody_from_obj(full_path, mesh.deltas, scaling);

        object->setUserPointer(&mesh);
    }
}

void scene::set_transducer(const btVector3 & position, const btVector3 & direction)
{
    transducer_pos = position;
    transducer_dir = direction;
}

template<unsigned int ray_count>
std::array<std::vector<ray_physics::segment>, ray_count> scene::cast_rays()
{
    using namespace ray_physics;

    std::array<std::vector<segment>, ray_count> segments;

    unsigned int tests = 0;
    unsigned int total_collisions = 0;

    ///step the simulation
    if (m_dynamicsWorld)
    {
        const float ray_start_step { 0.02f };

        for (auto & segments_vector : segments)
        {
            segments_vector.reserve(std::pow(2, ray::max_depth));
        }

        //#pragma omp parallel for
        for (size_t ray_i = 0; ray_i < ray_count; ray_i++)
        {
            auto & segments_vector = segments[ray_i];

            std::vector<ray> ray_stack;
            ray_stack.reserve(ray::max_depth-1);

            // Add first ray
            {
                ray first_ray
                {
                    transducer_pos + btVector3(0,0,ray_start_step * ray_i), // from
                    transducer_dir,                                         // initial direction
                    0,                                                      // depth
                    materials.at("GEL"),
                    nullptr,
                    initial_intensity,
                    frequency,
                    units::length::millimeter_t(0),                                                   // distance traveled
                    0                                                    // previous ray
                };
                ray_stack.push_back(first_ray);
            }

            while (ray_stack.size() > 0)
            {
                // Pop a ray from the stack and check if it collides

                auto & ray_ = ray_stack.at(ray_stack.size()-1);

                float r_length = ray_physics::max_ray_length(ray_);
                auto to = ray_.from + enlarge(ray_.direction, r_length);

                btCollisionWorld::ClosestRayResultCallback closestResults(ray_.from + 0.1f * ray_.direction,to);

                m_dynamicsWorld->rayTest(ray_.from + 0.1f * ray_.direction,to,closestResults);
                tests++;

                ray_stack.pop_back();

                if (closestResults.hasHit())
                {
                    // Substract ray intensity according to distance traveled
                    auto distance_before_hit = ray_.distance_traveled;
                    ray_physics::travel(ray_, distance_in_mm(ray_.from, closestResults.m_hitPointWorld));


                    if (ray_.depth < ray::max_depth)
                    {
                        // Calculate refraction and reflection directions and intensities

                        const auto organ = static_cast<mesh*>(closestResults.m_collisionObject->getUserPointer());
                        assert(organ);

                        auto result = ray_physics::hit_boundary(ray_, closestResults.m_hitPointWorld, closestResults.m_hitNormalWorld, *organ);

                        // Register collision creating a segment from the beggining of the ray to the collision point
                        segments_vector.emplace_back(segment{ray_.from, closestResults.m_hitPointWorld, ray_.direction, result.reflected_intensity, ray_.intensity, ray_.media.attenuation, distance_before_hit, ray_.media});

                        // Spawn reflection and refraction rays
                        if (result.refraction.intensity > ray::intensity_epsilon)
                        {
                            result.refraction.parent_collision = segments_vector.size()-1;
                            ray_stack.push_back(result.refraction);
                        }

                        if (result.reflection.intensity > ray::intensity_epsilon)
                        {
                            result.reflection.parent_collision = segments_vector.size()-1;
                            ray_stack.push_back(result.reflection);
                        }
                    }
                }
                else
                {
                    // Ray did not reach another media, add a data point at its end.
                    segments_vector.emplace_back(segment{ray_.from, to, ray_.direction, 0.0f, ray_.intensity, ray_.media.attenuation, ray_.distance_traveled + units::length::millimeter_t{r_length}, ray_.media});
                }
            }
        }

        for (auto & collision_vector : segments)
        {
            total_collisions += collision_vector.size();
        }
    }

    const float fps = 1.0f / (float( clock() - frame_start ) /  CLOCKS_PER_SEC);
    std::cout << fps << " " << tests << " " << total_collisions << std::endl;
    frame_start = clock();

    return segments;
}

void scene::parse_config(const nlohmann::json & config)
{
    working_dir = config.find("workingDirectory") != config.end() ? config.at("workingDirectory") : "";

    const auto & t_pos = config.at("transducerPosition");
    transducer_pos = {t_pos[0], t_pos[1], t_pos[2]};

    const auto & t_dir = config.at("transducerDirection");
    transducer_dir = {t_dir[0], t_dir[1], t_dir[2]};

    const auto & orig = config.at("origin");
    origin = {orig[0], orig[1], orig[2]};

    const auto & spac = config.at("spacing");
    spacing = {spac[0], spac[1], spac[2]};

    scaling = config.at("scaling");

    const auto & mats = config.at("materials");
    if (mats.is_array())
    {
        for (const auto & mat : mats)
        {
            materials[mat.at("name")] =
                {
                    mat.at("impedance"),
                    mat.at("attenuation"),
                    mat.at("mu0"),
                    mat.at("mu1"),
                    mat.at("sigma"),
                    mat.at("specularity")
                };
        }
    }
    else
    {
        throw std::runtime_error("materials must be an array");
    }

    const auto & meshes_ = config.at("meshes");
    if (meshes_.is_array())
    {
        for (const auto & mesh_ : meshes_)
        {
            const auto & deltas { mesh_.at("deltas") };
            meshes.emplace_back(mesh{
                        mesh_.at("file"),
                        mesh_.at("rigid"),
                        mesh_.at("vascular"),
                        {deltas[0], deltas[1], deltas[2]},
                        mesh_.at("outsideNormals"),
                        materials.at(mesh_.at("material")),
                        materials.at(mesh_.at("outsideMaterial"))});
        }
    }
    else
    {
        throw std::runtime_error("meshes must be an array");
    }
}

void scene::create_empty_world()
{
    m_collisionConfiguration = std::make_unique<btDefaultCollisionConfiguration>();

    m_dispatcher = std::make_unique<btCollisionDispatcher>(m_collisionConfiguration.get());

    m_broadphase = std::make_unique<btDbvtBroadphase>();

    m_solver = std::make_unique<btSequentialImpulseConstraintSolver>();

    m_dynamicsWorld = std::make_unique<btDiscreteDynamicsWorld>(m_dispatcher.get(),m_broadphase.get(),m_solver.get(),m_collisionConfiguration.get());

    m_dynamicsWorld->setGravity(btVector3(0,-10,0));
}

void scene::destroy_world()
{
    //delete collision shapes
    for (int j = 0; j < m_collisionShapes.size(); j++)
    {
        btCollisionShape* shape = m_collisionShapes[j];
        delete shape;
    }
    m_collisionShapes.clear();

    m_dynamicsWorld.reset();
    m_solver.reset();
    m_broadphase.reset();
    m_dispatcher.reset();
    m_collisionConfiguration.reset();
}

units::length::millimeter_t scene::distance_in_mm(const btVector3 & v1, const btVector3 & v2) const
{
    using namespace std;

    auto x_dist = abs(v1.getX() - v2.getX()) * spacing[0];
    auto y_dist = abs(v1.getY() - v2.getY()) * spacing[1];
    auto z_dist = abs(v1.getZ() - v2.getZ()) * spacing[2];

    return units::length::millimeter_t(sqrt(pow(x_dist,2) + pow(y_dist,2) + pow(z_dist,2)) * 10);
}

btVector3 scene::enlarge(const btVector3 & versor, float mm) const
{
    assert(versor.length2() < 1.1f);

    return mm/100.0f * btVector3 ( spacing[0] * versor.getX(),
                                   spacing[1] * versor.getY(),
                                   spacing[2] * versor.getZ() );
}

btRigidBody * scene::add_rigidbody_from_obj(const std::string & fileName, std::array<float, 3> deltas, float scaling)
{
    GLInstanceGraphicsShape* glmesh = load_mesh_from_obj(fileName, "");
    printf("[INFO] Obj loaded: Extracted %d verticed from obj file [%s]\n", glmesh->m_numvertices, fileName);

    const GLInstanceVertex& v = glmesh->m_vertices->at(0);
    btTriangleIndexVertexArray* tiva = new btTriangleIndexVertexArray(glmesh->m_numIndices / 3, &glmesh->m_indices->at(0), 3* sizeof(int),
                                                                      glmesh->m_numvertices, (btScalar*)(&(v.xyzw[0])), sizeof(GLInstanceVertex));

    btBvhTriangleMeshShape* shape = new btBvhTriangleMeshShape(tiva, true);

    m_collisionShapes.push_back(shape);

    float _scaling[4] = {scaling,scaling,scaling,1};

    btVector3 localScaling(_scaling[0],_scaling[1],_scaling[2]);
    shape->setLocalScaling(localScaling);

    btTransform startTransform;
    startTransform.setIdentity();

    //std::array<float, 3> origin { -18, -22, -5 }; // origin for organs scene
    float pos[4] = {deltas[0]*_scaling[0]*_scaling[0],deltas[1]*_scaling[1]*_scaling[1],deltas[2]*_scaling[2]*_scaling[2],0};
    btVector3 position(pos[0] + origin[0], pos[1] + origin[1], pos[2] + origin[2]);
    startTransform.setOrigin(position);

    btScalar mass(0.f);
    btVector3 localInertia(0, 0, 0);
    auto myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);
    auto body = new btRigidBody(cInfo);
    body->setUserIndex(-1);
    m_dynamicsWorld->addRigidBody(body);
    return body;
}

void scene::step(float delta_time)
{
    m_dynamicsWorld->stepSimulation(delta_time);
}

// TODO: Is this equals to distance_in_mm?
units::length::millimeter_t scene::distance(const btVector3 & from, const btVector3 & to) const
{
    // TODO: Use scaling in this calculation
    return units::length::millimeter_t(from.distance(to)*10.0f);
}

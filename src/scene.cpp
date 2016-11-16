#include "scene.h"

#include "objloader.h"
#include "ray.h"

#include <cmath>
#include <iostream>
#include <cassert>

template std::array<std::vector<ray_physics::segment>, 256> scene::cast_rays<256>();

scene::scene()
{
    create_empty_world();
}

scene::~scene()
{
    destroy_world();
}

void scene::init()
{
    struct mesh
    {
        std::string filename;
        bool is_rigid;
        std::array<float,3> deltas;
        material_id material;
    };

//    std::array<mesh,11> meshes
//    {{
//        {"aorta.obj", true, {152.533512115f, 174.472991943f, 105.106495678f}, material_id::BLOOD},
//        {"bones.obj", true, {188.265544891f, 202.440551758f, 105.599998474f}, material_id::BONE},
//        {"liver.obj", true, {141.238292694f, 176.429901123f, 130.10585022f}, material_id::LIVER},
//        {"cava.obj", true,  {206.332504272f, 192.29649353f, 104.897496045f}, material_id::BLOOD},
//        {"right_kidney.obj", true, {118.23374939f, 218.907501221f, 53.6022927761f}, material_id::KIDNEY},
//        {"left_kidney.obj", true,  {251.052993774f, 227.63949585f, 64.8468027115f}, material_id::KIDNEY},
//        {"right_suprarrenal.obj", true, {152.25050354f, 213.971496582f, 115.338005066f}, material_id::SUPRARRENAL},
//        {"left_suprarrenal.obj", true,  {217.128997803f, 209.525497437f, 102.477149963f}, material_id::SUPRARRENAL},
//        {"gallbladder.obj", true, {128.70715332f, 146.592498779f, 112.361503601f}, material_id::GALLBLADDER},
//        {"skin.obj", true,  {188.597551346f, 199.367202759f, 105.622316509f}, material_id::BONE},
//        {"porta.obj", true, {182.364089966f, 177.214996338f, 93.0034988523f}, material_id::BLOOD}
//    }};

    std::array<mesh,2> meshes
    {{
        {"BOX.obj", true, {0,0,0}, material_id::LIVER},
        {"SPHERE.obj", true, {0,0,0}, material_id::BONE}
    }};

    for (const auto & mesh : meshes)
    {
        const auto full_path = working_dir + mesh.filename;

        auto object = add_rigidbody_from_obj(full_path, mesh.deltas, 1.0f);

        // TODO: check the lifetime of this thing
        auto properties = new organ_properties(mesh.material);
        object->setUserPointer(properties);
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
        btVector3 initial_pos(-14,1.2,-3);
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
                    initial_pos + btVector3(0,0,ray_start_step * ray_i), // from
                    {1, 0, 0},                                           // initial direction
                    0,                                                   // depth
                    materials[material_id::GEL],
                    initial_intensity,
                    frequency,
                    0,                                                   // distance traveled
                    0                                                    // previous ray
                };
                ray_stack.push_back(first_ray);
            }

            while (ray_stack.size() > 0)
            {
                // Pop a ray from the stack and check if it collides

                auto ray_ = ray_stack.at(ray_stack.size()-1);

                float r_length = ray_physics::max_ray_length(ray_);
                auto to = ray_.from * 1.02f + enlarge(ray_.direction, r_length);

                btCollisionWorld::ClosestRayResultCallback closestResults(ray_.from,to);

                m_dynamicsWorld->rayTest(ray_.from,to,closestResults);
                tests++;

                ray_stack.pop_back();

                if (closestResults.hasHit())
                {
                    organ_properties * organ = static_cast<organ_properties*>(closestResults.m_collisionObject->getUserPointer());
                    const auto & organ_material = organ ? materials[organ->mat_id] : ray_.media;


                    // Substract ray intensity according to distance traveled
                    auto distance_before_hit = ray_.distance_traveled;
                    ray_physics::travel(ray_, distance_in_mm(ray_.from, closestResults.m_hitPointWorld));


                    if (ray_.depth < ray::max_depth)
                    {
                        // Calculate refraction and reflection directions and intensities

                        auto result = ray_physics::hit_boundary(ray_, closestResults.m_hitPointWorld, closestResults.m_hitNormalWorld, organ_material);

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
                    //segments_vector.emplace_back(segment{ray_.from, to, ray_.direction, 0.0f, ray_.intensity, ray_.media.attenuation, ray_.distance_traveled + r_length, ray_.media);
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

float scene::distance_in_mm(const btVector3 & v1, const btVector3 & v2) const
{
    using namespace std;

    auto x_dist = abs(v1.getX() - v2.getX()) * spacing[0];
    auto y_dist = abs(v1.getY() - v2.getY()) * spacing[1];
    auto z_dist = abs(v1.getZ() - v2.getZ()) * spacing[2];

    return sqrt(pow(x_dist,2) + pow(y_dist,2) + pow(z_dist,2));
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
    std::array<float, 3> origin { 0, 0, 0 };
    float pos[4] = {deltas[0]*_scaling[0]*_scaling[0],deltas[1]*_scaling[1]*_scaling[1],deltas[2]*_scaling[2]*_scaling[2],0};
    btVector3 position(pos[0] + origin[0], pos[1] + origin[1], pos[2] + origin[2]);
    startTransform.setOrigin(position);

    btScalar	mass(0.f);
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

float scene::distance(const btVector3 & from, const btVector3 & to)
{
    // TODO: Use scaling in this calculation
    return static_cast<float>(from.distance(to));
}

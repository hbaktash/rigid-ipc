#include "collision_constraint.hpp"

namespace ccd {
// !important: this needs to be define in the enum namespace

namespace opt {

    CollisionConstraint::CollisionConstraint(const std::string& name)
        : detection_method(HASH_GRID)
        , trajectory_type(SCREWING)
        , name_(name)
    {
    }

    void CollisionConstraint::settings(const nlohmann::json& json)
    {
        detection_method = json["detection_method"].get<DetectionMethod>();
        trajectory_type = json["trajectory_type"].get<TrajectoryType>();
    }

    nlohmann::json CollisionConstraint::settings() const
    {
        nlohmann::json json;
        json["detection_method"] = detection_method;
        json["trajectory_type"] = trajectory_type;
        return json;
    }

    void CollisionConstraint::construct_collision_set(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double> poses,
        const physics::Poses<double> displacements,
        ConcurrentImpacts& impacts) const
    {
        ccd::detect_collisions(
            bodies, poses, displacements, dim_to_collision_type(bodies.dim()),
            impacts, detection_method, trajectory_type);
    }

} // namespace opt
} // namespace ccd

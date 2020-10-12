#include "distance_barrier_constraint.hpp"

#include <igl/slice_mask.h>

#include <ipc/ipc.hpp>

#include <ccd/rigid_body_collision_detection.hpp>
#include <ccd/time_of_impact.hpp>
#include <geometry/distance.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    NLOHMANN_JSON_SERIALIZE_ENUM(
        BarrierType,
        { { BarrierType::IPC, "ipc" },
          { BarrierType::POLY_LOG, "poly_log" },
          { BarrierType::SPLINE, "spline" } })

    DistanceBarrierConstraint::DistanceBarrierConstraint(
        const std::string& name)
        : CollisionConstraint(name)
        , initial_barrier_activation_distance(1e-2)
        , barrier_type(BarrierType::IPC)
        , m_barrier_activation_distance(0.0)
    {
    }

    void DistanceBarrierConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        initial_barrier_activation_distance =
            json["initial_barrier_activation_distance"].get<double>();
        barrier_type = json["barrier_type"].get<BarrierType>();
    }

    nlohmann::json DistanceBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["initial_barrier_activation_distance"] =
            initial_barrier_activation_distance;
        json["barrier_type"] = barrier_type;
        return json;
    }

    void DistanceBarrierConstraint::initialize()
    {
        m_barrier_activation_distance = initial_barrier_activation_distance;
        CollisionConstraint::initialize();
    }

    bool DistanceBarrierConstraint::has_active_collisions(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1) const
    {
        PROFILE_POINT("collisions_detection");
        NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

        PROFILE_START();
        // This function will profile itself
        Candidates candidates;
        detect_collision_candidates(
            bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
            candidates, detection_method, trajectory_type,
            /*inflation_radius=*/0);

        PROFILE_START(NARROW_PHASE)
        bool has_collisions = has_active_collisions_narrow_phase(
            bodies, poses_t0, poses_t1, candidates);
        PROFILE_END(NARROW_PHASE)
        PROFILE_END();

        return has_collisions;
    }

    bool DistanceBarrierConstraint::has_active_collisions_narrow_phase(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const Candidates& candidates) const
    {
        for (const auto& ev_candidate : candidates.ev_candidates) {
            double toi, alpha;
            bool are_colliding = detect_edge_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, ev_candidate, toi, alpha,
                trajectory_type);
            if (are_colliding) {
                return true;
            }
        }
        for (const auto& fv_candidate : candidates.fv_candidates) {
            double toi, u, v;
            bool are_colliding = detect_face_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, fv_candidate, toi, u, v,
                trajectory_type);
            if (are_colliding) {
                return true;
            }
        }
        for (const auto& ee_candidate : candidates.ee_candidates) {
            double toi, edge0_alpha, edge1_alpha;
            bool are_colliding = detect_edge_edge_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, ee_candidate, toi, edge0_alpha,
                edge1_alpha, trajectory_type);
            if (are_colliding) {
                return true;
            }
        }
        return false;
    }

    double DistanceBarrierConstraint::compute_earliest_toi(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1) const
    {
        PROFILE_POINT("collisions_detection");
        NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

        PROFILE_START();
        // This function will profile itself
        Candidates candidates;
        detect_collision_candidates(
            bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
            candidates, detection_method, trajectory_type,
            /*inflation_radius=*/0);

        PROFILE_START(NARROW_PHASE)
        double earliest_toi = compute_earliest_toi_narrow_phase(
            bodies, poses_t0, poses_t1, candidates);
        PROFILE_END(NARROW_PHASE)
        PROFILE_END();

        return earliest_toi;
    }

    double DistanceBarrierConstraint::compute_earliest_toi_narrow_phase(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const Candidates& candidates) const
    {
        int collision_count = 0;
        double earliest_toi = 1;

        for (const auto& ev_candidate : candidates.ev_candidates) {
            double toi = std::numeric_limits<double>::infinity(), alpha;
            bool are_colliding = detect_edge_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, ev_candidate, toi, alpha,
                trajectory_type, earliest_toi);
            if (are_colliding && toi < earliest_toi) {
                collision_count++;
                earliest_toi = toi;
            }
        }
        for (const auto& fv_candidate : candidates.fv_candidates) {
            double toi = std::numeric_limits<double>::infinity(), u, v;
            bool are_colliding = detect_face_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, fv_candidate, toi, u, v,
                trajectory_type, earliest_toi);
            if (are_colliding && toi < earliest_toi) {
                collision_count++;
                earliest_toi = toi;
            }
        }
        for (const auto& ee_candidate : candidates.ee_candidates) {
            double toi = std::numeric_limits<double>::infinity(), edge0_alpha,
                   edge1_alpha;
            bool are_colliding = detect_edge_edge_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, ee_candidate, toi, edge0_alpha,
                edge1_alpha, trajectory_type, earliest_toi);
            if (are_colliding && toi < earliest_toi) {
                collision_count++;
                earliest_toi = toi;
            }
        }

        spdlog::info(
            "num_candidates={:d} num_collisions={:d} percentage={:g}%",
            candidates.size(), collision_count,
            candidates.size() == 0
                ? 0
                : double(collision_count) / candidates.size() * 100);

        return collision_count ? earliest_toi
                               : std::numeric_limits<double>::infinity();
    }

    void DistanceBarrierConstraint::construct_constraint_set(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        ipc::Constraints& constraint_set) const
    {
        // TODO: Update brute force version
        assert(detection_method == DetectionMethod::HASH_GRID);

        PROFILE_POINT("distance_barrier__construct_constraint_set");
        PROFILE_START();

        Eigen::MatrixXd V = bodies.world_vertices(poses);
        ipc::construct_constraint_set(
            /*V_rest=*/V, V, bodies.m_edges, bodies.m_faces,
            /*dhat=*/m_barrier_activation_distance, constraint_set,
            /*ignore_internal_vertices=*/false,
            /*vertex_group_ids=*/bodies.group_ids());

        PROFILE_END();
    }

    double DistanceBarrierConstraint::compute_minimum_distance(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses) const
    {
        ipc::Constraints constraint_set;
        construct_constraint_set(bodies, poses, constraint_set);
        Eigen::MatrixXd V = bodies.world_vertices(poses);
        return sqrt(ipc::compute_minimum_distance(
            V, bodies.m_edges, bodies.m_faces, constraint_set));
    }

} // namespace opt
} // namespace ccd

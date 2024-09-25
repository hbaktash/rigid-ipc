#include "distance_barrier_constraint.hpp"

#include <mutex>
#include <tbb/parallel_for_each.h>

#include <igl/slice_mask.h>
#include <ipc/ipc.hpp>

#include <ccd/rigid/broad_phase.hpp>
#include <ccd/rigid/rigid_body_hash_grid.hpp>
#include <ccd/save_queries.hpp>
#include <geometry/distance.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>
#include <profiler.hpp>

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"
#include "libqhull/qhull_a.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using orgQhull::Qhull;
using orgQhull::QhullFacet;
using orgQhull::QhullVertex;
using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace ipc::rigid {


std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, Eigen::MatrixXd>
get_convex_hull(Eigen::MatrixXd point_set){
    const size_t num_points = point_set.rows();
    const size_t dim = 3;
    char flags[64];
    // sprintf(flags, "qhull Qt");
    std::vector<double> pset_flat_vec(num_points*dim);
    for (size_t i = 0; i < num_points ; i++){
        pset_flat_vec[3*i]     = point_set(i, 0);
        pset_flat_vec[3*i + 1] = point_set(i, 1);
        pset_flat_vec[3*i + 2] = point_set(i, 2);
    }
    coordT* data = new coordT[num_points * dim];
    std::copy(pset_flat_vec.begin(), pset_flat_vec.end(), data);
    Qhull qhull("n", 3, num_points, data, "QJ Q3");
    
    // make data structures
    size_t my_count = 0, max_ind = 0;
    std::vector<size_t> indices;
    for (QhullVertex qvertex: qhull.vertexList()){
        my_count++;
        indices.push_back(qvertex.id()-1);
        if (qvertex.id() - 1 > max_ind)
            max_ind = qvertex.id() - 1;
    }
    if (qhull.vertexCount() > max_ind + 1)
        printf(" vert count is %d, my count: %d, max_ind: %d \n", qhull.vertexCount(), my_count, max_ind);
    std::vector<std::vector<size_t>> hull_faces;
    Eigen::MatrixXd poses(qhull.vertexCount(),3);  // 
    std::sort(indices.begin(), indices.end());
    std::map<size_t, size_t> re_indexing;
    for (size_t i = 0; i < indices.size(); i++){
        re_indexing[indices[i]] = i;
    }

    // tracing old vertices
    std::vector<size_t> old_indices(qhull.vertexCount());
    for (QhullVertex qvertex: qhull.vertexList()){
        size_t old_ind = qvertex.point().id();
        old_indices[re_indexing[qvertex.id() - 1]] = old_ind;
        // printf("old index %d new index- %d\n", old_ind, re_indexing[qvertex.id() - 1]);
    }
    for (QhullFacet qfacet: qhull.facetList().toStdVector()){
        std::vector<size_t> hull_face;
        // std::cout << qfacet << "\n";
        for (QhullVertex qvertex: qfacet.vertices().toStdVector()){
            // size_t id = qvertex.id() - 1; // qhull is 1-based for some reason
            size_t id = re_indexing[qvertex.id() - 1]; // qhull is 1-based for some reason
            hull_face.push_back(id);
            
            // position
            double x, y, z;
            const realT *c= qvertex.point().coordinates();
            for(int k= qvertex.dimension(); k--; ){
                // std::cout << " "<< k << " " << ; // QH11010 FIX: %5.2g
                if (k == 2) x = *c++;
                if (k == 1) y = *c++;
                if (k == 0) z = *c++;
            }
            poses.row(id) << x, y, z;
        }
        hull_faces.push_back(hull_face);
    }
    
    // build and return mesh
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    SurfaceMesh* surf_mesh = new SurfaceMesh(hull_faces);
    surf_mesh->greedilyOrientFaces();

    // make sure of outward orientation; can't just make qhull do it ??!
    Face f0 = surf_mesh->face(0);
    Vertex v0 = f0.halfedge().tailVertex(),
           v1 = f0.halfedge().tipVertex(),
           v2 = f0.halfedge().next().tipVertex();
    Eigen::Vector3d p0 = poses.row(v0.getIndex()),
                    p1 = poses.row(v1.getIndex()),
                    p2 = poses.row(v2.getIndex());
    Eigen::Vector3d p_other = poses.row(f0.halfedge().twin().next().tipVertex().getIndex()); // damn you qhull!!!!
    if ((p_other - p0).dot(cross(p1 - p0, p2 - p1)) > 0.){
        surf_mesh->invertOrientation(f0);
        surf_mesh->greedilyOrientFaces(); // will start with f0!
    }
    surf_mesh->compress();

    return {surf_mesh->getFaceVertexList(), old_indices, poses};
}


Eigen::VectorXi get_hull_indicator(Eigen::MatrixXd V, Eigen::MatrixXi F){
    int nf = F.rows(),
        nv = V.rows();
    Eigen::MatrixXd object_V = V.block(0, 0, nv - 4, 3);
    // get hull
    
    std::vector<std::vector<size_t>> hull_faces;
    std::vector<size_t> hull_to_input_map;
    Eigen::MatrixXd hull_poses;
    std::tie(hull_faces, hull_to_input_map, hull_poses) = get_convex_hull(object_V); // TODO
    Eigen::VectorXi is_on_hull(V.rows());
    is_on_hull.setZero();

    // index navigator
    for (std::vector<size_t> face: hull_faces){ 
        for (size_t hull_ind: face){
            int old_ind = hull_to_input_map[hull_ind]; // since ground indices come after object indices, we dont need a secon mapping
            is_on_hull(old_ind) = 1;
        }
    }
    Eigen::MatrixXi hull_faces_eigen(hull_faces.size(), 3);
    for (size_t i = 0; i < hull_faces.size(); i++){
        hull_faces_eigen.row(i) << hull_faces[i][0], hull_faces[i][1], hull_faces[i][2];
    }
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(hull_poses, hull_faces_eigen);
    viewer.data().set_face_based(true);
    viewer.launch();
    
    return is_on_hull;
}


NLOHMANN_JSON_SERIALIZE_ENUM(
    BarrierType,
    { { BarrierType::IPC, "ipc" },
      { BarrierType::POLY_LOG, "poly_log" },
      { BarrierType::SPLINE, "spline" } })

DistanceBarrierConstraint::DistanceBarrierConstraint(const std::string& name)
    : CollisionConstraint(name)
    , initial_barrier_activation_distance(1e-3)
    , barrier_type(BarrierType::IPC)
    , minimum_separation_distance(0.0)
    , m_barrier_activation_distance(0.0)
{
}

void DistanceBarrierConstraint::settings(const nlohmann::json& json)
{
    CollisionConstraint::settings(json);
    initial_barrier_activation_distance =
        json["initial_barrier_activation_distance"];
    minimum_separation_distance = json["minimum_separation_distance"];
    barrier_type = json["barrier_type"];
}

nlohmann::json DistanceBarrierConstraint::settings() const
{
    nlohmann::json json = CollisionConstraint::settings();
    json["initial_barrier_activation_distance"] =
        initial_barrier_activation_distance;
    json["minimum_separation_distance"] = minimum_separation_distance;
    json["barrier_type"] = barrier_type;
    return json;
}

void DistanceBarrierConstraint::initialize()
{
    m_barrier_activation_distance = initial_barrier_activation_distance;
    CollisionConstraint::initialize();
}

bool DistanceBarrierConstraint::has_active_collisions(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1) const
{
    PROFILE_POINT("DistanceBarrierConstraint::has_active_collisions");
    NAMED_PROFILE_POINT(
        "DistanceBarrierConstraint::has_active_collisions_narrow_phase",
        NARROW_PHASE);

    PROFILE_START();
    // This function will profile itself
    Candidates candidates;
    detect_collision_candidates(
        bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
        candidates, detection_method, trajectory_type,
        /*inflation_radius=*/minimum_separation_distance / 2.0);

    PROFILE_START(NARROW_PHASE)
    bool has_collisions = has_active_collisions_narrow_phase(
        bodies, poses_t0, poses_t1, candidates);
    PROFILE_END(NARROW_PHASE)
    PROFILE_END();

    return has_collisions;
}

bool DistanceBarrierConstraint::has_active_collisions_narrow_phase(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const Candidates& candidates) const
{
    TrajectoryType overloaded_trajectory =
        trajectory_type == TrajectoryType::PIECEWISE_LINEAR
        ? TrajectoryType::RIGID
        : trajectory_type;

    for (const auto& ev_candidate : candidates.ev_candidates) {
        double toi;
        bool are_colliding = edge_vertex_ccd(
            bodies, poses_t0, poses_t1, ev_candidate, toi,
            overloaded_trajectory);
        if (are_colliding) {
            // save_ccd_candidate(bodies, poses_t0, poses_t1, ev_candidate);
            return true;
        }
    }
    for (const auto& fv_candidate : candidates.fv_candidates) {
        double toi;
        bool are_colliding = face_vertex_ccd(
            bodies, poses_t0, poses_t1, fv_candidate, toi,
            overloaded_trajectory);
        if (are_colliding) {
            // save_ccd_candidate(bodies, poses_t0, poses_t1, fv_candidate);
            return true;
        }
    }
    for (const auto& ee_candidate : candidates.ee_candidates) {
        double toi;
        bool are_colliding = edge_edge_ccd(
            bodies, poses_t0, poses_t1, ee_candidate, toi,
            overloaded_trajectory);
        if (are_colliding) {
            // save_ccd_candidate(bodies, poses_t0, poses_t1, ee_candidate);
            return true;
        }
    }
    return false;
}

double DistanceBarrierConstraint::compute_earliest_toi(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1) const
{
    PROFILE_POINT("DistanceBarrierConstraint::compute_earliest_toi");
    PROFILE_START();
    // This function will profile itself
    Candidates candidates;
    detect_collision_candidates(
        bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
        candidates, detection_method, trajectory_type,
        /*inflation_radius=*/minimum_separation_distance / 2.0);

    double earliest_toi = compute_earliest_toi_narrow_phase(
        bodies, poses_t0, poses_t1, candidates);
    PROFILE_END();

    return earliest_toi;
}

double DistanceBarrierConstraint::compute_earliest_toi_narrow_phase(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const Candidates& candidates) const
{
    NAMED_PROFILE_POINT(
        "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase",
        NARROW_PHASE);
    // NOTE: These are disabled because they are not thread safe
    // NAMED_PROFILE_POINT(
    //     "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase:"
    //     "edge_edge",
    //     FV_NARROW_PHASE);
    // NAMED_PROFILE_POINT(
    //     "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase:"
    //     "face_vertex",
    //     EE_NARROW_PHASE);
    // NAMED_PROFILE_POINT(
    //     "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase:"
    //     "edge_vertex",
    //     EV_NARROW_PHASE);

    PROFILE_START(NARROW_PHASE);

    int collision_count = 0;
    double earliest_toi = 1;
    std::mutex earliest_toi_mutex;

    const size_t num_ev = candidates.ev_candidates.size();
    const size_t num_ee = candidates.ee_candidates.size();
    const size_t num_fv = candidates.fv_candidates.size();

    // Do a single block range over all three candidate vectors
    tbb::parallel_for(
        tbb::blocked_range<int>(0, candidates.size()),
        [&](tbb::blocked_range<int> r) {
            for (int i = r.begin(); i < r.end(); i++) {
                double toi = std::numeric_limits<double>::infinity();
                bool are_colliding;

                if (i < num_ev) {
                    // PROFILE_START(EV_NARROW_PHASE);
                    are_colliding = edge_vertex_ccd(
                        bodies, poses_t0, poses_t1, candidates.ev_candidates[i],
                        toi, trajectory_type, earliest_toi,
                        minimum_separation_distance);
                    // PROFILE_END(EV_NARROW_PHASE);
                } else if (i - num_ev < num_ee) {
                    // PROFILE_START(EE_NARROW_PHASE);
                    are_colliding = edge_edge_ccd(
                        bodies, poses_t0, poses_t1,
                        candidates.ee_candidates[i - num_ev], toi,
                        trajectory_type, earliest_toi,
                        minimum_separation_distance);
                    // PROFILE_END(EE_NARROW_PHASE);
                } else {
                    assert(i - num_ev - num_ee < num_fv);
                    // PROFILE_START(FV_NARROW_PHASE);
                    are_colliding = face_vertex_ccd(
                        bodies, poses_t0, poses_t1,
                        candidates.fv_candidates[i - num_ev - num_ee], toi,
                        trajectory_type, earliest_toi,
                        minimum_separation_distance);
                    // PROFILE_END(FV_NARROW_PHASE);
                }

                if (are_colliding && toi == 0) {
                    if (i < num_ev) {
                        spdlog::error("Edge-vertex CCD resulted in toi=0!");
                        save_ccd_candidate(
                            bodies, poses_t0, poses_t1,
                            candidates.ev_candidates[i]);
                    } else if (i - num_ev < num_ee) {
                        spdlog::error("Edge-edge CCD resulted in toi=0!");
                        save_ccd_candidate(
                            bodies, poses_t0, poses_t1,
                            candidates.ee_candidates[i - num_ev]);
                    } else {
                        assert(i - num_ev - num_ee < num_fv);
                        spdlog::error("Face-vertex CCD resulted in toi=0!");
                        save_ccd_candidate(
                            bodies, poses_t0, poses_t1,
                            candidates.fv_candidates[i - num_ev - num_ee]);
                    }
                }

                if (are_colliding) {
                    std::scoped_lock lock(earliest_toi_mutex);
                    collision_count++;
                    if (toi < earliest_toi) {
                        earliest_toi = toi;
                    }
                }
            }
        });

    double percent_correct = candidates.size() == 0
        ? 100
        : (double(collision_count) / candidates.size() * 100);
    PROFILE_MESSAGE(
        NARROW_PHASE, "num_candidates,num_collisions,percentage",
        fmt::format(
            "{:d},{:d},{:g}%", candidates.size(), collision_count,
            percent_correct));

    spdlog::debug(
        "num_candidates={:d} num_collisions={:d} percentage={:g}%",
        candidates.size(), collision_count, percent_correct);

    PROFILE_END(NARROW_PHASE);

    return collision_count ? earliest_toi
                           : std::numeric_limits<double>::infinity();
}

void DistanceBarrierConstraint::construct_constraint_set(
    const RigidBodyAssembler& bodies,
    const PosesD& poses,
    Constraints& constraint_set) const
{
    static PosesD cached_poses;
    static Constraints cached_constraint_set;

    if (bodies.num_bodies() <= 1) {
        return;
    }

    if (poses == cached_poses) {
        constraint_set = cached_constraint_set;
        return;
    }

    PROFILE_POINT("DistanceBarrierConstraint::construct_constraint_set");
    PROFILE_START();

    const double& dhat = m_barrier_activation_distance;
    const double& dmin = minimum_separation_distance;
    const double inflation_radius = (dhat + dmin) / 2.0;

    Candidates candidates;
    detect_collision_candidates_rigid(
        bodies, poses, dim_to_collision_type(bodies.dim()), candidates,
        detection_method, inflation_radius);
    
    // TODO remove from candids

    std:: cout << " here @@@@@@@@@@@@@@@@@@@@@ ------------------- @@@@@@@@@@@@@@@@@@@@@ " << std::endl;
    
    Eigen::MatrixXd V1 = bodies.world_vertices(poses);
    Eigen::MatrixXi F = bodies.m_faces;
    Eigen::VectorXi is_on_hull = get_hull_indicator(V1, F);
    
    
    // candidates.ev_candidates.clear();
    // candidates.ee_candidates.clear();
    // for (FaceVertexCandidate cand: candidates.fv_candidates){
    //     cand.face_index;
    // }


    Eigen::MatrixXd V = bodies.world_vertices(poses);
    ipc::construct_constraint_set(
        candidates, /*V_rest=*/V, V, bodies.m_edges, bodies.m_faces,
        /*dhat=*/dhat, constraint_set, bodies.m_faces_to_edges,
        /*dmin=*/dmin);

    PROFILE_END();

    cached_poses = poses;
    cached_constraint_set = constraint_set;
}

double DistanceBarrierConstraint::compute_minimum_distance(
    const RigidBodyAssembler& bodies, const PosesD& poses) const
{
    PROFILE_POINT("DistanceBarrierConstraint::compute_minimum_distance");
    PROFILE_START();

    Constraints constraint_set;
    construct_constraint_set(bodies, poses, constraint_set);
    Eigen::MatrixXd V = bodies.world_vertices(poses);
    double minimum_distance = sqrt(ipc::compute_minimum_distance(
        V, bodies.m_edges, bodies.m_faces, constraint_set));

    PROFILE_END();

    return minimum_distance;
}

} // namespace ipc::rigid

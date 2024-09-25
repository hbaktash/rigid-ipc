#include "rigid_body.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <autodiff/autodiff_types.hpp>
#include <finitediff.hpp>
#include <logger.hpp>
#include <physics/mass.hpp>
#include <profiler.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/flatten.hpp>
#include <utils/not_implemented_error.hpp>

// #include <igl/Timer.h>
#include <igl/edges.h>

#include <Eigen/Core>
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"
#include "libqhull/qhull_a.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using orgQhull::Qhull;
using orgQhull::QhullFacet;
using orgQhull::QhullVertex;


namespace ipc::rigid {

std::tuple<Eigen::MatrixXi, std::vector<size_t>, Eigen::MatrixXd>
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
    std::vector<std::vector<size_t>> hull_faces;
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
    //
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

    hull_faces = surf_mesh->getFaceVertexList(); // oriented
    
    Eigen::MatrixXi F(hull_faces.size(), 3);
    for (size_t i = 0; i < hull_faces.size(); i++){
        F.row(i) << hull_faces[i][0], hull_faces[i][1], hull_faces[i][2];
    }

    return {F, old_indices, poses};
}


// Eigen::VectorXi get_hull_indicator(Eigen::MatrixXd V){
//     int nv = V.rows();
//     Eigen::MatrixXd object_V = V.block(0, 0, nv - 4, 3);
//     // get hull
    
//     std::vector<std::vector<size_t>> hull_faces;
//     std::vector<size_t> hull_to_input_map;
//     Eigen::MatrixXd hull_poses;
//     std::tie(hull_faces, hull_to_input_map, hull_poses) = get_convex_hull(object_V); // TODO
//     Eigen::VectorXi is_on_hull(V.rows());
//     is_on_hull.setZero();

//     // index navigator
//     for (std::vector<size_t> face: hull_faces){ 
//         for (size_t hull_ind: face){
//             int old_ind = hull_to_input_map[hull_ind]; // since ground indices come after object indices, we dont need a secon mapping
//             is_on_hull(old_ind) = 1;
//         }
//     }
//     Eigen::MatrixXi hull_faces_eigen(hull_faces.size(), 3);
//     for (size_t i = 0; i < hull_faces.size(); i++){
//         hull_faces_eigen.row(i) << hull_faces[i][0], hull_faces[i][1], hull_faces[i][2];
//     }
//     // igl::opengl::glfw::Viewer viewer;
//     // viewer.data().set_mesh(hull_poses, hull_faces_eigen);
//     // viewer.data().set_face_based(true);
//     // viewer.launch();
    
//     return is_on_hull;
// }

void center_vertices(
    Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    PoseD& pose)
{
    int dim = vertices.cols();

    vertices.rowwise() += pose.position.transpose();

    // compute the center of mass several times to get more accurate
    pose.position.setZero(dim);
    for (int i = 0; i < 10; i++) {
        double mass;
        VectorMax3d com;
        MatrixMax3d inertia;
        compute_mass_properties(
            vertices, dim == 2 || faces.size() == 0 ? edges : faces, mass, com,
            inertia);
        vertices.rowwise() -= com.transpose();
        pose.position += com;
        if (com.squaredNorm() < 1e-8) {
            break;
        }
    }
}

RigidBody::RigidBody(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const PoseD& pose,
    const PoseD& velocity,
    const PoseD& force,
    const double density,
    const VectorMax6b& is_dof_fixed,
    const bool oriented,
    const int group_id,
    const RigidBodyType type,
    const double kinematic_max_time,
    const std::deque<PoseD>& kinematic_poses)
    : group_id(group_id)
    , type(type)
    , vertices(vertices)
    , edges(edges)
    , faces(faces)
    , is_dof_fixed(is_dof_fixed)
    , is_oriented(oriented)
    , mesh_selector(vertices.rows(), edges, faces)
    , pose(pose)
    , velocity(velocity)
    , force(force)
    , kinematic_max_time(kinematic_max_time)
    , kinematic_poses(kinematic_poses)
{
    assert(dim() == pose.dim());
    assert(dim() == velocity.dim());
    assert(dim() == force.dim());
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    if (type == RigidBodyType::STATIC) {
        this->is_dof_fixed.setOnes(this->is_dof_fixed.size());
    } else if (this->is_dof_fixed.array().all()) {
        this->type = RigidBodyType::STATIC;
    }

    center_vertices(this->vertices, edges, faces, this->pose);
    VectorMax3d center_of_mass;
    MatrixMax3d I;
    compute_mass_properties(
        this->vertices,
        dim() == 2 || faces.size() == 0 ? edges : faces, //
        mass, center_of_mass, I);
    // assert(center_of_mass.squaredNorm() < 1e-8);

    // Mass above is actually volume in m³ and density is Kg/m³
    mass *= density;
    if (dim() == 3) {
        // Got this from Chrono: https://bit.ly/2RpbTl1
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
        double threshold = I.lpNorm<Eigen::Infinity>() * 1e-16;
        I = (threshold < I.array().abs()).select(I, 0.0);
        es.compute(I);
        if (es.info() != Eigen::Success) {
            spdlog::error("Eigen decompostion of the inertia tensor failed!");
        }
        moment_of_inertia = density * es.eigenvalues();
        if ((moment_of_inertia.array() < 0).any()) {
            spdlog::warn(
                "Negative moment of inertia ({}), inverting.",
                fmt_eigen(moment_of_inertia));
            // Avoid negative epsilon inertias
            moment_of_inertia =
                (moment_of_inertia.array() < 0)
                    .select(-moment_of_inertia, moment_of_inertia);
        }
        R0 = es.eigenvectors();
        // Ensure that we have an orientation preserving transform
        if (R0.determinant() < 0.0) {
            R0.col(0) *= -1.0;
        }
        assert(R0.isUnitary(1e-9));
        assert(fabs(R0.determinant() - 1.0) <= 1.0e-9);
        int num_rot_dof_fixed =
            is_dof_fixed.tail(PoseD::dim_to_rot_ndof(dim())).count();
        if (num_rot_dof_fixed == 2) {
            // Convert moment of inertia to world coordinates
            // https://physics.stackexchange.com/a/268812
            moment_of_inertia = -I.diagonal().array() + I.diagonal().sum();
            R0.setIdentity();
        } else if (num_rot_dof_fixed == 1) {
            spdlog::warn("Rigid body dynamics with two rotational DoF has "
                         "not been tested thoroughly.");
        }
        // R = RᵢR₀
        Eigen::AngleAxisd r = Eigen::AngleAxisd(
            Eigen::Matrix3d(this->pose.construct_rotation_matrix() * R0));
        this->pose.rotation = r.angle() * r.axis();
        // v = Rv₀ + p = RᵢR₀v₀ + p = RᵢR₀R₀ᵀv₀ + p
        this->vertices = this->vertices * R0; // R₀ᵀ * V₀ᵀ = V₀ * R₀
        // ω = R₀ᵀω₀ (ω₀ expressed in body coordinates)
        this->velocity.rotation = R0.transpose() * this->velocity.rotation;
        Eigen::Matrix3d Q_t0 = this->pose.construct_rotation_matrix();
        this->Qdot = Q_t0 * Hat(this->velocity.rotation);
        // τ = R₀ᵀτ₀ (τ₀ expressed in body coordinates)
        // NOTE: this transformation will be done later
        // this->force.rotation = R0.transpose() * this->force.rotation;


        // Hull stuff
        // only modifying V, F, E
        //  ; other geometric properties are computed after this
        //  ; center of mass is already computed
        if (this->vertices.rows() > 4){ // not ground
            std::cout<< " ######### modifying object to its convex hull ######### "<< std::endl;
            std::vector<size_t> hull_to_input_map;
            Eigen::MatrixXi hull_F;
            Eigen::MatrixXd hull_V;
            std::tie(hull_F, hull_to_input_map, hull_V) = get_convex_hull(this->vertices); // TODO
            this->vertices = hull_V;
            this->faces = hull_F;
            std::cout<< " original vertices: "<< vertices.rows() << " original faces: "<< faces.rows() << std::endl;
            std::cout<< " hull vertices: "<< hull_V.rows() << " hull faces: "<< hull_F.rows() << std::endl;
            // just following read_obj.cpp; not sure why conservative resize is used.
            Eigen::MatrixXi faceE, tmp_E;
            igl::edges(this->faces, faceE); 
            tmp_E.conservativeResize(tmp_E.rows() + faceE.rows(), 2);
            tmp_E.bottomRows(faceE.rows()) = faceE;
            this->edges = tmp_E; 
            std::cout<<" original edges: "<< edges.rows() << " \n hull edges: "<< tmp_E.rows() << std::endl;
            mesh_selector = MeshSelector(this->vertices.rows(), this->edges, this->faces);
        }
        else {
            std::cout<< " @@@@@@@@@@ ground probably @@@@@@@@@@ "<< std::endl;
        }
    } else {
        moment_of_inertia = density * I.diagonal();
        R0 = Eigen::Matrix<double, 1, 1>::Identity();
    }

    // Zero out the velocity and forces of fixed dof
    this->velocity.zero_dof(is_dof_fixed, R0);
    this->force.zero_dof(is_dof_fixed, R0);

    // Update the previous pose and velocity to reflect the changes made
    // here
    this->pose_prev = this->pose;
    this->velocity_prev = this->velocity;

    this->acceleration = PoseD::Zero(dim());
    this->Qddot.setZero();

    // Compute and construct some useful constants
    mass_matrix.resize(ndof());
    mass_matrix.diagonal().head(pos_ndof()).setConstant(mass);
    mass_matrix.diagonal().tail(rot_ndof()) = moment_of_inertia;

    r_max = this->vertices.rowwise().norm().maxCoeff();

    average_edge_length = 0;
    for (long i = 0; i < this->edges.rows(); i++) {
        average_edge_length +=
            (this->vertices.row(this->edges(i, 0)) - this->vertices.row(this->edges(i, 1)))
                .norm();
    }
    if (this->edges.rows() > 0) {
        average_edge_length /= this->edges.rows();
    }
    assert(std::isfinite(average_edge_length));

    init_bvh();
}

void RigidBody::init_bvh()
{
    PROFILE_POINT("RigidBody::init_bvh");
    PROFILE_START();

    // heterogenous bounding boxes
    std::vector<std::array<Eigen::Vector3d, 2>> aabbs(
        num_codim_vertices() + num_codim_edges() + num_faces());

    for (size_t i = 0; i < num_codim_vertices(); i++) {
        size_t vi = mesh_selector.codim_vertices_to_vertices(i);
        if (dim() == 2) {
            aabbs[i][0][2] = 0;
            aabbs[i][1][2] = 0;
        }
        aabbs[i][0].head(dim()) = vertices.row(i);
        aabbs[i][1].head(dim()) = vertices.row(i);
    }

    size_t start_i = num_codim_vertices();
    for (size_t i = 0; i < num_codim_edges(); i++) {
        size_t ei = mesh_selector.codim_edges_to_edges(i);
        const auto& e0 = vertices.row(edges(ei, 0));
        const auto& e1 = vertices.row(edges(ei, 1));

        if (dim() == 2) {
            aabbs[start_i + i][0][2] = 0;
            aabbs[start_i + i][1][2] = 0;
        }
        aabbs[start_i + i][0].head(dim()) = e0.cwiseMin(e1);
        aabbs[start_i + i][1].head(dim()) = e0.cwiseMax(e1);
    }

    start_i += num_codim_edges();
    for (size_t i = 0; i < num_faces(); i++) {
        assert(dim() == 3);
        const auto& f0 = vertices.row(faces(i, 0));
        const auto& f1 = vertices.row(faces(i, 1));
        const auto& f2 = vertices.row(faces(i, 2));
        aabbs[start_i + i][0] = f0.cwiseMin(f1).cwiseMin(f2);
        aabbs[start_i + i][1] = f0.cwiseMax(f1).cwiseMax(f2);
    }

    bvh.init(aabbs);

    PROFILE_END();
}

Eigen::MatrixXd RigidBody::world_velocities() const
{
    // compute ẋ = Q̇ * x_B + q̇
    // where Q̇ = Qω̂
    if (dim() == 2) {
        MatrixMax3d Q_dt =
            pose.construct_rotation_matrix() * Hat(velocity.rotation);
        return (vertices * Q_dt.transpose()).rowwise()
            + velocity.position.transpose();
    }
    return (vertices * Qdot.transpose()).rowwise()
        + velocity.position.transpose();
}

void RigidBody::compute_bounding_box(
    const PoseD& pose_t0,
    const PoseD& pose_t1,
    VectorMax3d& box_min,
    VectorMax3d& box_max) const
{
    PROFILE_POINT("RigidBody::compute_bounding_box");
    PROFILE_START();

    // If the body is not rotating then just use the linearized
    // trajectory
    if (type == RigidBodyType::STATIC
        || (pose_t0.rotation.array() == pose_t1.rotation.array()).all()) {
        Eigen::MatrixXd V0 = world_vertices(pose_t0);
        box_min = V0.colwise().minCoeff();
        box_max = V0.colwise().maxCoeff();
        Eigen::MatrixXd V1 = world_vertices(pose_t1);
        box_min = box_min.cwiseMin(V1.colwise().minCoeff().transpose());
        box_max = box_max.cwiseMax(V1.colwise().maxCoeff().transpose());
    } else {
        // Use the maximum radius of the body to bound all rotations
        box_min = pose_t0.position.cwiseMin(pose_t1.position).array() - r_max;
        box_max = pose_t0.position.cwiseMax(pose_t1.position).array() + r_max;
    }

    PROFILE_END();
}

} // namespace ipc::rigid

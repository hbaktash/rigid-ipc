#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <autodiff/autodiff.h>
#include <physics/rigid_body.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace physics {

    class RigidBodyAssembler {
    public:
        RigidBodyAssembler() {}
        ~RigidBodyAssembler() {}

        /// @brief inits assembler to use this set of rigid-bodies
        void init(const std::vector<RigidBody>& rbs);

        // World Vertices Functions
        // --------------------------------------------------------------------

        template <typename T>
        Eigen::MatrixX<T> world_vertices(
            const std::vector<Eigen::MatrixXX3<T>>& rotations,
            const std::vector<Eigen::VectorX3<T>>& positions) const;
        template <typename T>
        Eigen::MatrixX<T> world_vertices(const Poses<T>& poses) const;
        template <typename T>
        Eigen::MatrixX<T> world_vertices(const Eigen::VectorX<T>& dofs) const
        {
            return world_vertices(Pose<T>::dofs_to_poses(dofs, dim()));
        }

        Eigen::MatrixXd world_vertices(
            const RigidBody::Step step = RigidBody::CURRENT_STEP) const;
        Eigen::MatrixXd world_vertices_t0() const
        {
            return world_vertices(RigidBody::PREVIOUS_STEP);
        }
        Eigen::MatrixXd world_vertices_t1() const
        {
            return world_vertices(RigidBody::CURRENT_STEP);
        }

        template <typename T>
        Eigen::VectorX3<T>
        world_vertex(const Pose<T>& pose, const int vertex_idx) const;
        template <typename T>
        Eigen::VectorX3<T>
        world_vertex(const Poses<T>& poses, const int vertex_idx) const;

        Eigen::MatrixXd world_velocities() const;

        void world_vertices_gradient(
            const Poses<double>& poses,
            Eigen::SparseMatrix<double>& grad) const;
        void world_vertices_gradient(
            const Eigen::VectorXd& dof, Eigen::SparseMatrix<double>& grad) const
        {
            return world_vertices_gradient(
                Pose<double>::dofs_to_poses(dof, dim()), grad);
        }

        void global_to_local_vertex(
            const long global_vertex_id,
            long& rigid_body_id,
            long& local_vertex_id) const;
        void global_to_local_edge(
            const long global_edge_id,
            long& rigid_body_id,
            long& local_edge_id) const;
        void global_to_local_face(
            const long global_face_id,
            long& rigid_body_id,
            long& local_face_id) const;

        // Ridig Body CM Functions
        // --------------------------------------------------------------------
        /// @brief assemble rigid body poses to a single vector
        Poses<double> rb_poses(const bool previous = false) const;
        Poses<double> rb_poses_t0() const { return rb_poses(true); }
        Poses<double> rb_poses_t1() const { return rb_poses(false); }
        /// @brief set rigid body poses
        void set_rb_poses(const Poses<double>& poses);

        // --------------------------------------------------------------------

        long num_vertices() const { return m_body_vertex_id.back(); }
        long num_edges() const { return m_body_edge_id.back(); }
        long num_faces() const { return m_body_face_id.back(); }
        size_t num_bodies() const { return m_rbs.size(); }
        int dim() const { return m_rbs.size() ? m_rbs[0].dim() : 0; }

        long vertex_id_to_body_id(long vi) const
        {
            return m_vertex_to_body_map(vi);
        }
        long edge_id_to_body_id(long ei) const
        {
            return m_vertex_to_body_map(m_edges(ei, 0));
        }
        long face_id_to_body_id(long fi) const
        {
            return m_vertex_to_body_map(m_faces(fi, 0));
        }

        const RigidBody& vertex_id_to_body(long vi) const
        {
            return m_rbs[vertex_id_to_body_id(vi)];
        }
        const RigidBody& edge_id_to_body(long ei) const
        {
            return m_rbs[edge_id_to_body_id(ei)];
        }
        const RigidBody& face_id_to_body(long fi) const
        {
            return m_rbs[face_id_to_body_id(fi)];
        }

        const Eigen::VectorXi& group_ids() const { return m_vertex_group_ids; }

        /// Get a vector of body ids where each body is close to at least one
        /// other body.
        std::vector<int> close_bodies(
            const physics::Poses<double>& poses_t0,
            const physics::Poses<double>& poses_t1,
            const double inflation_radius) const;
        std::vector<int> close_bodies_brute_force(
            const physics::Poses<double>& poses_t0,
            const physics::Poses<double>& poses_t1,
            const double inflation_radius) const;

        std::vector<RigidBody> m_rbs;

        /// @brief indexes for body
        std::vector<long> m_body_vertex_id;
        std::vector<long> m_body_face_id;
        std::vector<long> m_body_edge_id;

        double average_edge_length; ///< Average edge length

        /// @brief average value of mass matrix diagonal (including inertia)
        double average_mass;

        /// @brief indexes for vertices
        Eigen::VectorXi m_vertex_to_body_map;

        /// @brief re-indexed edges of the whole system
        Eigen::MatrixXi m_edges;

        /// @brief re-indexed faces of the whole system
        Eigen::MatrixXi m_faces;

        /// @brief mass_matrix of the rigid bodies dof
        Eigen::DiagonalMatrixXd m_rb_mass_matrix;

        /// @brief flags for rb degrees of freedom
        Eigen::VectorXb is_rb_dof_fixed;

        /// @brief flag for vertices degrees of freedom (used for visualization)
        Eigen::MatrixXb is_dof_fixed;

    protected:
        /// @brief Group ids per vertex
        Eigen::VectorXi m_vertex_group_ids;
    };

} // namespace physics
} // namespace ccd

#include "rigid_body_assembler.tpp"

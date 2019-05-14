// Methods for optimizing the displacments with a non-linear interference volume
// constraint.
#include <opt/displacement_opt.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/not_implemented_error.hpp>

#include <autodiff/finitediff.hpp>

#include <logger.hpp>

namespace ccd {
namespace opt {

    ParticlesDisplProblem::ParticlesDisplProblem()
        : constraint(nullptr)
    {
    }

    ParticlesDisplProblem::~ParticlesDisplProblem() {}

    void ParticlesDisplProblem::initialize(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2i& E, const Eigen::MatrixX2d& U,
        CollisionConstraint& cstr)
    {
        vertices = V;
        edges = E;
        displacements = U;

        constraint = &cstr;
        constraint->initialize(V, E, U);
        initProblem();
    }

    void ParticlesDisplProblem::initProblem()
    {
        u_ = displacements;
        u_.resize(displacements.size(), 1);

        num_vars = int(u_.size());
        x0.resize(num_vars);
        x_lower.resize(num_vars);
        x_upper.resize(num_vars);
        x_lower.setConstant(NO_LOWER_BOUND);
        x_upper.setConstant(NO_UPPER_BOUND);

        // TODO: this is wrong, it will depend on constraint type
        num_constraints = int(edges.rows());
        g_lower.resize(num_constraints);
        g_upper.resize(num_constraints);
        g_lower.setConstant(0.0);
        g_upper.setConstant(NO_UPPER_BOUND);
    }

    double ParticlesDisplProblem::eval_f(const Eigen::VectorXd& x)
    {
        return (x - u_).squaredNorm() / 2.0;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        return (x - u_);
    }

    Eigen::MatrixXd ParticlesDisplProblem::eval_hessian_f(
        const Eigen::VectorXd& x)
    {
        return Eigen::MatrixXd::Identity(x.size(), x.size());
    }

    Eigen::SparseMatrix<double> ParticlesDisplProblem::eval_hessian_f_sparse(
        const Eigen::VectorXd& x)
    {
        Eigen::SparseMatrix<double> A(int(x.size()), int(x.size()));
        A.setIdentity();
        return A;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_g(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        Eigen::VectorXd gx;
        if (constraint->recompute_collision_set) {
            constraint->detecteCollisions(vertices, edges, Uk);
        }
        constraint->compute_constraints(vertices, edges, Uk, gx);
        return gx;
    };

    Eigen::MatrixXd ParticlesDisplProblem::eval_jac_g(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        Eigen::MatrixXd jac_gx;
        if (constraint->recompute_collision_set) {
            constraint->detecteCollisions(vertices, edges, Uk);
        }
        constraint->compute_constraints_jacobian(vertices, edges, Uk, jac_gx);

        return jac_gx;
    };

    std::vector<Eigen::MatrixXd> ParticlesDisplProblem::eval_hessian_g(
        const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        std::vector<Eigen::MatrixXd> hess_gx;
        if (constraint->recompute_collision_set) {
            constraint->detecteCollisions(vertices, edges, Uk);
        }
        constraint->compute_constraints_hessian(vertices, edges, Uk, hess_gx);
        return hess_gx;
    };


} // namespace opt
} // namespace ccd

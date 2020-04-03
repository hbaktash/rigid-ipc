#include "pose.hpp"

#include <Eigen/Geometry>
#include <tbb/tbb.h>

#include <logger.hpp>
#include <utils/is_zero.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace physics {

    template <typename T>
    Pose<T>::Pose()
        : position()
        , rotation()
    {
    }

    template <typename T>
    Pose<T>::Pose(
        const Eigen::VectorX3<T>& position, const Eigen::VectorX3<T>& rotation)
        : position(position)
        , rotation(rotation)
    {
    }

    template <typename T> Pose<T>::Pose(const Eigen::VectorX6<T>& dof)
    {
        if (dof.size() == dim_to_ndof(2)) {
            position = dof.head(dim_to_pos_ndof(2));
            rotation = dof.tail(dim_to_rot_ndof(2));
        } else if (dof.size() == dim_to_ndof(3)) {
            position = dof.head(dim_to_pos_ndof(3));
            rotation = dof.tail(dim_to_rot_ndof(3));
        } else {
            throw NotImplementedError("Unknown pose convertion for given ndof");
        }
    }

    template <typename T>
    Pose<T>::Pose(const T& x, const T& y, const T& theta)
        : Pose(Eigen::Vector2<T>(x, y), Eigen::Vector1<T>())
    {
        rotation << theta;
    }

    template <typename T>
    Pose<T>::Pose(
        const T& x,
        const T& y,
        const T& z,
        const T& theta_x,
        const T& theta_y,
        const T& theta_z)
        : Pose(
              Eigen::Vector3<T>(x, y, z),
              Eigen::Vector3<T>(theta_x, theta_y, theta_z))
    {
    }

    template <typename T> Pose<T> Pose<T>::Zero(int dim)
    {
        assert(dim == 2 || dim == 3);
        return Pose(
            Eigen::VectorX<T>::Zero(Pose<T>::dim_to_pos_ndof(dim)),
            Eigen::VectorX<T>::Zero(Pose<T>::dim_to_rot_ndof(dim)));
    }

    template <typename T>
    Poses<T> Pose<T>::dofs_to_poses(const Eigen::VectorX<T>& dofs, int dim)
    {
        int ndof = dim_to_ndof(dim);
        int num_poses = dofs.size() / ndof;
        assert(dofs.size() % ndof == 0);
        Poses<T> poses;
        poses.reserve(num_poses);
        for (int i = 0; i < num_poses; i++) {
            poses.emplace_back(dofs.segment(i * ndof, ndof));
        }
        return poses;
    }

    template <typename T>
    Eigen::VectorX<T> Pose<T>::poses_to_dofs(const Poses<T>& poses)
    {
        const int ndof = poses.size() ? poses[0].ndof() : 0;
        Eigen::VectorX<T> dofs(poses.size() * ndof);
        for (size_t i = 0; i < poses.size(); i++) {
            assert(poses[i].ndof() == ndof);
            dofs.segment(i * ndof, ndof) = poses[i].dof();
        }
        return dofs;
    }

    template <typename T> Eigen::VectorX6<T> Pose<T>::dof() const
    {
        Eigen::VectorX6<T> pose_dof(ndof());
        pose_dof.head(pos_ndof()) = position;
        pose_dof.tail(rot_ndof()) = rotation;
        return pose_dof;
    }

    // Replace a selected dof with the dof in other.
    template <typename T>
    void Pose<T>::select_dof(
        const Eigen::VectorX6b& is_dof_selected,
        const Pose<T>& other,
        const Eigen::MatrixXX3d& R)
    {
        assert(is_dof_selected.size() == this->ndof());
        assert(other.dim() == this->dim());
        // R should be a rotation
        assert(R.isUnitary(1e-9));
        assert(fabs(R.determinant() - 1.0) < 1.0e-6);
        position =
            is_dof_selected.head(pos_ndof()).select(other.position, position);
        rotation = R.transpose()
            * is_dof_selected.tail(rot_ndof())
                  .select(R * other.rotation, R * rotation);
    }

    // Zero out the i-th dof if is_dof_zero(i) == true.
    template <typename T>
    void Pose<T>::zero_dof(
        const Eigen::VectorX6b& is_dof_zero, const Eigen::MatrixXX3d& R)
    {
        select_dof(is_dof_zero, Pose<T>::Zero(dim()), R);
    }

    template <typename T>
    Eigen::MatrixXX3<T> Pose<T>::construct_rotation_matrix() const
    {
        if (dim() == 2) {
            return Eigen::Rotation2D<T>(rotation(0)).toRotationMatrix();
        } else {
            T angle = rotation.norm();
            Eigen::Vector3<T> axis =
                rotation / (is_zero(angle) ? T(1.0) : angle);
            return Eigen::AngleAxis<T>(angle, axis).toRotationMatrix();
        }
    }

    template <typename T>
    Eigen::Quaternion<T> Pose<T>::construct_quaternion() const
    {
        if (dim() == 2) {
            return Eigen::Quaternion<T>(
                Eigen::AngleAxis<T>(rotation(0), Eigen::Vector3<T>::UnitZ()));
        } else {
            T angle = rotation.norm();
            Eigen::Vector3<T> axis = rotation;
            if (is_zero(angle)) {
                // The error of doing this will be on the order of the
                // adjustment.
                angle += T(1e-16);
            }
            axis /= angle;
            return Eigen::Quaternion<T>(Eigen::AngleAxis<T>(angle, axis));
        }
    }

    template <typename T>
    Pose<T>
    Pose<T>::interpolate(const Pose<T>& pose0, const Pose<T>& pose1, T t)
    {
        static_assert(
            !std::is_same<T, Interval>::value,
            "You cannot use interpolate with Intervals");
        assert(pose0.dim() == pose1.dim());
        Pose<T> pose;
        pose.position = (pose1.position - pose0.position) * t + pose0.position;
        if (pose0.dim() == 2) {
            pose.rotation =
                (pose1.rotation - pose0.rotation) * t + pose0.rotation;
        } else {
            Eigen::Quaternion<T> q0 = pose0.construct_quaternion();
            Eigen::Quaternion<T> q1 = pose1.construct_quaternion();
            Eigen::Quaternion<T> q = q0.slerp(t, q1);
            Eigen::AngleAxis<T> aa(q);
            pose.rotation = aa.angle() * aa.axis();
        }
        return pose;
    }

    template <typename T> bool Pose<T>::operator==(const Pose<T>& other) const
    {
        return this->position == other.position
            && this->rotation == other.rotation;
    }

    template <typename T> Pose<T>& Pose<T>::operator*=(const T& x)
    {
        this->position *= x;
        this->rotation *= x;
        return *this;
    }

    template <typename T> Pose<T> Pose<T>::operator/(const T& x) const
    {
        return Pose<T>(this->position / x, this->rotation / x);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Operations on vector of Poses

    template <typename T>
    Poses<T> interpolate(const Poses<T>& poses0, const Poses<T>& poses1, T t)
    {
        Poses<T> poses(poses0.size());
        tbb::parallel_for(size_t(0), poses.size(), [&](size_t i) {
            poses[i] = Pose<T>::interpolate(poses0[i], poses1[i], t);
        });
        return poses;
    }

    template <typename T> Poses<T> operator*(const Poses<T>& poses, const T& x)
    {
        Poses<T> product = poses;
        tbb::parallel_for(
            size_t(0), product.size(), [&](size_t i) { product[i] *= x; });
        return product;
    }

    template <typename T1, typename T2>
    Poses<T2> cast(const Poses<T1>& poses_T1)
    {
        Poses<T2> poses_T2;
        poses_T2.reserve(poses_T1.size());
        for (int i = 0; i < poses_T1.size(); i++) {
            poses_T2.push_back(poses_T1[i].template cast<T2>());
        }
        return poses_T2;
    }

} // namespace physics

template <typename T> Eigen::Matrix3<T> rotate_to_z(Eigen::Vector3<T> n)
{
    if (n.norm() == T(0)) {
        return Eigen::Matrix3<T>::Identity();
    }
    return Eigen::Quaternion<T>::FromTwoVectors(n, Eigen::Vector3<T>::UnitZ())
        .toRotationMatrix();
}

template <typename T> Eigen::Matrix3<T> rotate_around_z(const T& theta)
{
    Eigen::Matrix3<T> R;
    R.row(0) << cos(theta), -sin(theta), T(0);
    R.row(1) << sin(theta), cos(theta), T(0);
    R.row(2) << T(0), T(0), T(1);
    return R;
}

template <typename T>
void decompose_to_z_screwing(
    const physics::Pose<T>& pose_t0,
    const physics::Pose<T>& pose_t1,
    Eigen::Matrix3<T>& R0,
    Eigen::Matrix3<T>& P,
    T& omega)
{
    // Decompose the inbetween rotation as a rotation around the z-axis:
    //     R = Pᵀ R_z P
    // Where R = R₁R₀ᵀ, P is a rotation from n̂ to ẑ, and R_z is a rotation
    // of θ around the z-axis.
    R0 = pose_t0.construct_rotation_matrix();
    Eigen::Matrix3<T> R1 = pose_t1.construct_rotation_matrix();
    Eigen::AngleAxis<T> r(R1 * R0.transpose());
    omega = r.angle();
    P = rotate_to_z(r.axis());
}

} // namespace ccd

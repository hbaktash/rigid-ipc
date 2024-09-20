/************************************************************************
*
* ADOBE CONFIDENTIAL
* ___________________
*
* Copyright [first year code created] Adobe
* All Rights Reserved.
*
* NOTICE: All information contained herein is, and remains
* the property of Adobe and its suppliers, if any. The intellectual
* and technical concepts contained herein are proprietary to Adobe
* and its suppliers and are protected by all applicable intellectual
* property laws, including trade secret and copyright laws.
* Dissemination of this information or reproduction of this material
* is strictly forbidden unless prior written permission is obtained
* from Adobe.
*************************************************************************
*/
#include "unsupported/Eigen/EulerAngles"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "nlohmann/json.hpp"

#include <SimState.hpp>
#include <igl/opengl/glfw/Viewer.h>
// #include "polyscope/nlohmann/json.hpp"
// #include "bullet3/examples/BasicExample.h"
// #include "args.hxx"

// #include "ipc/ipc.hpp"

using namespace geometrycentral;
using namespace geometrycentral::surface;


// GC stuff
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;


int main(int argc, char* argv[])
{
    // Load a mesh
    std::tie(mesh_ptr, geometry_ptr) = readManifoldSurfaceMesh("../meshes/spot.obj");
    mesh = mesh_ptr.release();
    geometry = geometry_ptr.release();

    ipc::rigid::SimState sim;
    bool success = sim.load_scene("../fixtures/3D/examples/example.json");
    if (!success) {
        return 1;
    }
    
    // std::string fout = "../fixtures/3D/examples/my_wrapper_sim.json";
    // // sim.run_simulation(fout);

    // run steps manually
    // spdlog::info("Starting simulation ", scene_file);
    // spdlog::info("Running {} iterations", m_max_simulation_steps);

    igl::Timer timer;

    Eigen::MatrixXd V1 = sim.problem_ptr.get()->vertices();
    Eigen::MatrixXi F1 = sim.problem_ptr.get()->faces();
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V1, F1);
    viewer.data().set_face_based(true);
    // Draw dfdU as little white lines from each vertex
    viewer.launch();

    timer.start();
    sim.m_solve_collisions = true;
    // print_progress_bar(0, m_max_simulation_steps, 0);
    for (int i = 0; i < sim.m_max_simulation_steps; ++i) {
        sim.simulation_step();
        sim.save_simulation_step();
        // spdlog::info("Finished it={} sim_step={}", i + 1, m_num_simulation_steps);
    }
    std::cout << "simulation ended with steps: "<< sim.state_sequence.size() << std::endl;
    timer.stop();
    std::cout << "that took " << timer.getElapsedTime() << " seconds" << std::endl;

    // Draw the mesh
    Eigen::MatrixXd V = sim.problem_ptr.get()->vertices();
    Eigen::MatrixXi F = sim.problem_ptr.get()->faces();
    // viewer.data().clear();
    igl::opengl::glfw::Viewer viewer2;
    viewer2.data().set_mesh(V, F);
    viewer2.data().set_face_based(true);
    // Draw dfdU as little white lines from each vertex
    viewer2.launch();

    return EXIT_SUCCESS;	
}
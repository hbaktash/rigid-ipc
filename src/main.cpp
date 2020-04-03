#include <igl/opengl/glfw/Viewer.h>

#include <logger.hpp>
#include <viewer/UISimState.hpp>

#include <io/json_to_mjcf.hpp>

int main(int argc, char* argv[])
{
    if (argc > 1) {
        int mode = std::stoi(argv[1]);
        switch (mode) {
        case 1: { // json to MJCF
            spdlog::info("transform json to MJCF");

            if (argc < 3) {
                spdlog::error("need input json file path!");
                return -1;
            }

            return json_to_mjcf(argv[2]);
        }

        case 0: // simulation
        default:
            break;
        }
    }

    // tbb::task_scheduler_init init(1);
    spdlog::set_level(spdlog::level::info);
    ccd::UISimState ui;
    ui.launch();
}

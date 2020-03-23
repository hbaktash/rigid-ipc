#include "time_stepper_factory.hpp"

#include <time_stepper/dmv_time_stepper.hpp>
#include <time_stepper/sympletic_euler_time_stepper.hpp>

namespace ccd {

const TimeStepperFactory& TimeStepperFactory::factory()
{
    static TimeStepperFactory instance;
    return instance;
}

TimeStepperFactory::TimeStepperFactory()
{
    using namespace time_stepper;
    time_steppers.emplace(
        SympleticEulerTimeStepper::default_name(),
        std::make_shared<SympleticEulerTimeStepper>());
    time_steppers.emplace(
        DMVTimeStepper::default_name(), std::make_shared<DMVTimeStepper>());
}

std::shared_ptr<time_stepper::TimeStepper>
TimeStepperFactory::get_time_stepper(const std::string& name) const
{
    auto it = time_steppers.find(name);
    assert(it != time_steppers.end());
    return it->second;
}

std::shared_ptr<time_stepper::TimeStepper>
TimeStepperFactory::get_default_time_stepper(int dim) const
{
    assert(dim == 2 || dim == 3);
    using namespace time_stepper;
    return dim == 2
        ? get_time_stepper(SympleticEulerTimeStepper::default_name())
        : get_time_stepper(DMVTimeStepper::default_name());
}

} // namespace ccd

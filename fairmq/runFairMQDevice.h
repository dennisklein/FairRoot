/********************************************************************************
 *    Copyright (C) 2017 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include <fairmq/DeviceRunner.h>
#include <boost/program_options.hpp>
#include <memory>
#include <string>

template <typename R>
class GenericFairMQDevice : public FairMQDevice
{
  public:
    GenericFairMQDevice(R func) : r(func) {}

  protected:
    virtual bool ConditionalRun() { return r(*static_cast<FairMQDevice*>(this)); }

  private:
    R r;
};

template <typename R>
FairMQDevice* makeDeviceWithConditionalRun(R r)
{
    return new GenericFairMQDevice<R>(r);
}

using FairMQDevicePtr = FairMQDevice*;

// to be implemented by the user to return a child class of FairMQDevice
FairMQDevicePtr getDevice(const FairMQProgOptions& config);

// to be implemented by the user to add custom command line options (or just with empty body)
void addCustomOptions(boost::program_options::options_description&);

int main(int argc, const char** argv)
{
    DeviceRunner runner;

    runner.Subscribe([](DeviceRunner::Hook hook, DeviceRunner& runner){
        switch (hook) {
            case Hook::LoadPlugins:
                // e.g. runner.fPluginManager->LoadPlugin("asdf");
                break;
            case Hook::SetCustomCmdLineOptions:
                boost::program_options::options_description customOptions("Custom options");
                addCustomOptions(customOptions);
                runner.fConfig.AddToCmdLineOptions(customOptions);
                break;
            case Hook::ModifyRawCmdLineArgs:
                // e.g. runner.fRawCmdLineArgs.push_back("--blubb");
                break;
            case Hook::InstantiateDevice:
                runner.fDevice = std::shared_ptr<FairMQDevice>{getDevice(runner.fConfig)};
                break;
            default:
                break;
        }
    });

    return runner.RunWithExceptionHandlers();

    // Run without builtin catch all exception handler, just:
    // return runner.Run();
}

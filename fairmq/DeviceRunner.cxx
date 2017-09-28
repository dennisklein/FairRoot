/********************************************************************************
 *    Copyright (C) 2017 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include "DeviceRunner.cxx"
#include <fairmq/Tools.h>

using namespace fair::mq;

DeviceRunner::DeviceRunner(int argc, char const* argv[])
: fRawCmdLineArgs{tools::ToStrVector(argc, argv)}
, fPluginManager{PluginManager::MakeFromCommandLineOptions(fRawCmdLineArgs)}
, fDevice{nullptr}
{
}

auto Subscribe(std::function<void(Step, DeviceRunner&)> hook) -> void
{
    fEvents.Subscribe<DeviceRunner&>("runner", hook);
}

auto Unsubscribe() -> void
{
    fEvents.Unsubscribe<DeviceRunner&>("runner");
}

auto DeviceRunner::Run() -> int
{
    fEvents.Emit<DeviceRunner&>(Hook::LoadPlugins, *this);

    // Load builtin plugins last
    pluginManager->LoadPlugin("s:control");

    fEvents.Emit<DeviceRunner&>(Hook::SetCustomCmdLineOptions, *this);

    pluginManager->ForEachPluginProgOptions([&fConfig](boost::program_options::options_description options){
        fConfig.AddToCmdLineOptions(options);
    });
    config.AddToCmdLineOptions(pluginManager->ProgramOptions());

    fEvents.Emit<DeviceRunner&>(Hook::ModifyCmdLineOptions, *this);

    // TODO implement vector interface
    fConfig.ParseAll(argc, argv, true);

    fEvents.Emit<DeviceRunner&>(Hook::InstantiateDevice, *this);

    if (!fDevice)
    {
        LOG(ERROR) << "getDevice(): no valid device provided. Exiting.";
        return 1;
    }

    // Handle --print-channels
    fDevice->RegisterChannelEndpoints();
    if (fConfig.Count("print-channels"))
    {
        fDevice->PrintRegisteredChannels();
        fDevice->ChangeState(FairMQDevice::END);
        return 0;
    }

    // Handle --version
    if (fConfig.Count("version"))
    {
        std::cout << "User device version: " << fDevice->GetVersion() << std::endl;
        std::cout << "FAIRMQ_INTERFACE_VERSION: " << FAIRMQ_INTERFACE_VERSION << std::endl;
        fDevice->ChangeState(FairMQDevice::END);
        return 0;
    }

    LOG(DEBUG) << "PID: " << getpid();

    // Configure device
    fDevice->SetConfig(fConfig);

    // Initialize plugin services
    fPluginManager->EmplacePluginServices(&config, device);

    // Instantiate and run plugins
    fPluginManager->InstantiatePlugins();

    // Wait for control plugin to release device control
    fPluginManager->WaitForPluginsToReleaseDeviceControl();

    return 0;
}

auto DeviceRunner::RunWithExceptionHandlers() -> int
{
    try
    {
        return Run();
    }
    catch (std::exception& e)
    {
        LOG(ERROR) << "Unhandled exception reached the top of main: " << e.what() << ", application will now exit";
        return 1;
    }
    catch (...)
    {
        LOG(ERROR) << "Non-exception instance being thrown. Please make sure you use std::runtime_exception() instead. Application will now exit.";
        return 1;
    }
}

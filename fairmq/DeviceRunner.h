/********************************************************************************
 *    Copyright (C) 2017 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#ifndef FAIR_MQ_DEVICERUNNER_H
#define FAIR_MQ_DEVICERUNNER_H

#include <fairmq/EventManager.h>
#include <fairmq/PluginManager.h>
#include <FairMQDevice.h>
#include <FairMQLogger.h>
#include <options/FairMQProgOptions.h>

#include <functional>
#include <string>
#include <vector>

namespace fair
{
namespace mq
{

/**
 * @class DeviceRunner DeviceRunner.h <fairmq/DeviceRunner.h>
 * @brief Utility class to facilitate a convenient top-level device launch/shutdown.
 *
 * Runs a FairMQ in the context with config and plugin support.
 *
 * For customization user hooks are executed at various steps during device launch/shutdown.
 */
class DeviceRunner
{
  public:
    enum class Step
    {
        LoadPlugins,
        SetCustomCmdLineOptions,
        ModifyRawCmdLineArgs,
        InstantiateDevice,
    };

    DeviceRunner(int argc, char const* argv[]);

    auto Run() -> int;
    auto RunWithExceptionHandlers() -> int;

    auto Subscribe(std::function<void(Step, DeviceRunner&)> hook) -> void;
    auto Unsubscribe() -> void;

    std::vector<std::string> fRawCmdLineArgs;
    std::unique_ptr<PluginManager> fPluginManager;
    FairMQProgOptions fConfig;
    std::shared_ptr<FairMQDevice> fDevice;

  private:
    struct UserHook : Event<Step> {};
    EventManager fEvents;
};

} /* namespace mq */
} /* namespace fair */

#endif /* FAIR_MQ_DEVICERUNNER_H */

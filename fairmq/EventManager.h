/********************************************************************************
 * Copyright (C) 2014-2017 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH  *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#ifndef FAIR_MQ_EVENTMANAGER_H
#define FAIR_MQ_EVENTMANAGER_H

#include <mutex>
#include <string>
#include <type_traits>
#include <utility>

#include <boost/signals2.hpp>

namespace fair
{
namespace mq
{

// Inherit from this base event type to create custom event types
template<typename K>
struct Event
{
    using KeyType = const K;
};

/**
 * @class EventManager EventManager.h <fairmq/EventManager.h>
 * @brief Manages event callbacks from different subscribers
 *
 * The event manager stores a set of callbacks and associates them with
 * events depending on the callback signature. The first callback
 * argument must be of a special key type determined by the event type.
 *
 * Callbacks can be subscribed/unsubscribed based on a subscriber id,
 * the event type, and the callback signature.
 *
 * Events can be emitted based on event type and callback signature.
 *
 * The event manager is thread-safe, as long as different threads use different subscriber ids or (inclusive OR) callback signatures.
 *
 * Double subscriptions are not supported.
 */
class EventManager
{
  public:
    template<typename E, typename ...Args>
    using Callback = std::function<void(typename E::KeyType, Args...)>;
    template<typename E, typename ...Args>
    using Signal = boost::signals2::signal<void(typename E::KeyType, Args...)>;

    template<typename E, typename ...Args>
    static auto Subscribe(const std::string& subscriber, Callback<E, Args...> callback) -> void
    {
        Scribe<0, E, Args...>(true, subscriber, callback);
    }

    template<typename E, typename ...Args>
    static auto Unsubscribe(const std::string& subscriber) -> void
    {
        Scribe<0, E, Args...>(false, subscriber);
    }

    template<typename E, typename ...Args>
    static auto Emit(typename E::KeyType& key, Args&&... args) -> void
    {
        GetSignal<E, Args...>()(key, std::forward<Args>(args)...);
    }

  private:
    static constexpr int fkMaxSubscriptions{100};

    // acts as a map from string -> (connection, E, Args) with linear search time
    template<int N, typename E, typename ...Args>
    static auto Scribe(bool Connecting, const std::string& subscriber, Callback<E, Args...> callback = {}) -> void
    {
        static std::string name{subscriber};

        if (name == subscriber)
        {
            static boost::signals2::connection connection;

            if (Connecting)
            {
                assert(!connection.connected());
                connection = GetSignal<E,Args...>().connect(callback);
            }
            else
            {
                assert(connection.connected());
                connection.disconnect();
            }
        }
        else
        {
            constexpr auto bounded = N > fkMaxSubscriptions ? fkMaxSubscriptions : N+1;
            assert(N <= fkMaxSubscriptions);

            Scribe<bounded, E, Args...>(Connecting, subscriber, callback);
        }
    }

    // acts as map from (E, Args) -> Signal
    template<typename E, typename ...Args>
    static auto GetSignal() -> Signal<E, Args...>&
    {
        static Signal<E, Args...> signal;

        return signal;
    }
}; /* class EventManager */

} /* namespace mq */
} /* namespace fair */

#endif /* FAIR_MQ_EVENTMANAGER_H */

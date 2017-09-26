/********************************************************************************
 *    Copyright (C) 2017 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include <gtest/gtest.h>
#include <fairmq/EventManager.h>
#include <string>

namespace
{

using namespace std;
using namespace fair::mq;

struct TestEvent : fair::mq::Event<std::string> {};

TEST(EventManager, Basics)
{
    int call_counter = 0;
    int call_counter2 = 0;
    int value = 0;
    string value2;

    EventManager::Callback<TestEvent, int> callback{
        [&](TestEvent::KeyType& key, int newValue){
            ++call_counter;
            if (key == "test") value = newValue;
        }
    };
    EventManager::Callback<TestEvent, string> callback2{
        [&](TestEvent::KeyType& key, string newValue){
            ++call_counter2;
            if (key == "test") value2 = newValue;
        }
    };

    EventManager::Subscribe<TestEvent, int>("foo_subscriber", callback);
    EventManager::Emit<TestEvent>(TestEvent::KeyType{"test"}, 42);
    ASSERT_EQ(call_counter, 1);
    ASSERT_EQ(value, 42);

    EventManager::Unsubscribe<TestEvent, int>("foo_subscriber");
    EventManager::Emit<TestEvent>(TestEvent::KeyType{"test"}, 13);
    ASSERT_EQ(call_counter, 1);

    EventManager::Subscribe<TestEvent, int>("foo_subscriber", callback);
    EventManager::Subscribe<TestEvent, string>("foo_subscriber", callback2);
    // two different callbacks:
    EventManager::Emit<TestEvent>(TestEvent::KeyType{"test"}, 9000);
    EventManager::Emit<TestEvent>(TestEvent::KeyType{"test"}, string{"over 9000"});
    ASSERT_EQ(call_counter, 2);
    ASSERT_EQ(value, 9000);
    ASSERT_EQ(call_counter2, 1);
    ASSERT_EQ(value2, "over 9000");
}

} // namespace

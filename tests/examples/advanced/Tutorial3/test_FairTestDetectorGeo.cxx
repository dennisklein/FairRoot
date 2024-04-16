/********************************************************************************
 *   Copyright (C) 2024 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH     *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include <FairTestDetectorGeo.h>
#include <catch2/catch_test_macros.hpp>
#include <string>

TEST_CASE("FairTestDetectorGeo")
{
    FairTestDetectorGeo geo;

    SECTION("Constructor")
    {
        // Check that members are correctly set by constructor
        std::string retVal = geo.GetName();
        REQUIRE(retVal == "torino");
        REQUIRE(geo.getMaxModules() == 10);
    }

    SECTION("getModuleName")
    {
        std::string retVal = geo.getModuleName(0);
        REQUIRE(retVal == "torino1");

        retVal = geo.getModuleName(1);
        REQUIRE(retVal == "torino2");

        retVal = geo.getModuleName(99);
        REQUIRE(retVal == "torino100");
    }

    SECTION("getEleName")
    {
        std::string retVal = geo.getEleName(0);
        REQUIRE(retVal == "torino1");

        retVal = geo.getEleName(1);
        REQUIRE(retVal == "torino2");

        retVal = geo.getEleName(99);
        REQUIRE(retVal == "torino100");
    }

    SECTION("getModNumInMod")
    {
        REQUIRE(geo.getModNumInMod("torino1") == 0);

        REQUIRE(geo.getModNumInMod("torino2") == 1);

        REQUIRE(geo.getModNumInMod("torino9") == 8);

        //      The current code doesn't work for numbers with more than one digit
        REQUIRE(geo.getModNumInMod("torino10") == 0);
    }
}

/********************************************************************************
 *   Copyright (C) 2024 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH     *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include <FairTutorialDet1Geo.h>
#include <catch2/catch_test_macros.hpp>

#include <string>

TEST_CASE("FairTutorialDet1Geo")
{
    FairTutorialDet1Geo geo;

    SECTION("Constructor")
    {
        // Check that members are correctly set by constructor
        std::string retVal = geo.GetName();
        REQUIRE(retVal == "tutdet");
        REQUIRE(geo.getMaxModules() == 10);
    }

    SECTION("getModuleName")
    {
        std::string retVal = geo.getModuleName(0);
        REQUIRE(retVal == "tutdet1");

        retVal = geo.getModuleName(1);
        REQUIRE(retVal == "tutdet2");

        retVal = geo.getModuleName(99);
        REQUIRE(retVal == "tutdet100");
    }

    SECTION("getEleName")
    {
        std::string retVal = geo.getEleName(0);
        REQUIRE(retVal == "tutdet1");

        retVal = geo.getEleName(1);
        REQUIRE(retVal == "tutdet2");

        retVal = geo.getEleName(99);
        REQUIRE(retVal == "tutdet100");
    }

    geo.getModNumInMod("tutdet1");

    SECTION("getModNumInMod")
    {
        REQUIRE(geo.getModNumInMod("tutdet1") == 0);

        REQUIRE(geo.getModNumInMod("tutdet2") == 1);

        REQUIRE(geo.getModNumInMod("tutdet9") == 8);

//      The current code doesn't work for numbers with more than one digit
        REQUIRE(geo.getModNumInMod("tutdet10") == 0);
    }

/*
    SECTION("test")
    {
        std::string retVal = geo.test("tutdet1");
        REQUIRE(retVal == "1");

        retVal = geo.test("tutdet2");
        REQUIRE(retVal == "2");

        retVal = geo.test("tutdet99");
        REQUIRE(retVal == "99");
    }
*/

}

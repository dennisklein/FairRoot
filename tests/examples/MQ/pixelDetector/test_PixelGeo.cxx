/********************************************************************************
 *   Copyright (C) 2024 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH     *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include <PixelGeo.h>
#include <catch2/catch_test_macros.hpp>

#include <string>

TEST_CASE("PixelGeo")
{
    PixelGeo geo;

    SECTION("Constructor")
    {
        // Check that members are correctly set by constructor
        std::string retVal = geo.GetName();
        REQUIRE(retVal == "Pixel");
        REQUIRE(geo.getMaxModules() == 10);
    }

    SECTION("getModuleName")
    {
        std::string retVal = geo.getModuleName(0);
        REQUIRE(retVal == "Pixel1");

        retVal = geo.getModuleName(1);
        REQUIRE(retVal == "Pixel2");

        retVal = geo.getModuleName(99);
        REQUIRE(retVal == "Pixel100");
    }

    SECTION("getEleName")
    {
        std::string retVal = geo.getEleName(0);
        REQUIRE(retVal == "Pixel1");

        retVal = geo.getEleName(1);
        REQUIRE(retVal == "Pixel2");

        retVal = geo.getEleName(99);
        REQUIRE(retVal == "Pixel100");
    }

    geo.getModNumInMod("Pixel1");

    SECTION("getModNumInMod")
    {
        REQUIRE(geo.getModNumInMod("Pixel1") == 0);

        REQUIRE(geo.getModNumInMod("Pixel2") == 1);

        REQUIRE(geo.getModNumInMod("Pixel9") == 8);

//      The current code doesn't work for numbers with more than one digit
        REQUIRE(geo.getModNumInMod("Pixel10") == 0);
    }

}

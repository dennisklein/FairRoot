// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DetectorsBase/Aligner.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include <chrono>

O2ParamImpl(o2::base::Aligner);

using namespace o2::base;

using DetID = o2::detectors::DetID;

DetID::mask_t Aligner::getDetectorsMask() const
{
  return DetID::getMask(mDetectors) & (~DetID::getMask(DetID::CTP));
}

void Aligner::applyAlignment(long timestamp, DetID::mask_t addMask) const
{
  DetID::mask_t msk = getDetectorsMask() & addMask;
  if (msk.none()) {
    return;
  }
  if (!gGeoManager) {
    throw std::runtime_error("Geometry is not loaded, cannot apply alignment");
  }
  if (gGeoManager->IsLocked()) {
    throw std::runtime_error("Geometry is locked, cannot apply alignment");
  }

  throw std::runtime_error("alignment from the CCDB is not part of this build");
}

long Aligner::getTimeStamp() const
{
  return mTimeStamp > 0 ? mTimeStamp : std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

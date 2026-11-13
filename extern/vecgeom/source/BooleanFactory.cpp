// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.

#include "VecGeom/management/BooleanFactory.h"
#include "VecGeom/volumes/UnplacedBooleanVolume.h"
#include "VecGeom/volumes/UnplacedMultiUnion.h"
#include "VecGeom/volumes/UnplacedScaledShape.h"

#include <vector>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
namespace {

struct Component {
  VUnplacedVolume const *solid;
  Transformation3D transform;
};

// Only associative unions can be expanded into this component list. In
// particular, a subtraction on the right remains a single cutter expression.
void CollectUnion(VUnplacedVolume const *solid, Transformation3D const &frame, std::vector<Component> &components)
{
  auto const *boolean = BooleanHelper::GetBooleanStruct(solid);
  if (boolean && boolean->fOp == kUnion) {
    for (auto const *placed : {boolean->fLeftVolume, boolean->fRightVolume}) {
      Transformation3D transform(frame);
      transform.MultiplyFromRight(*placed->GetTransformation());
      CollectUnion(placed->GetUnplacedVolume(), transform, components);
    }
  } else if (auto const *multi = dynamic_cast<UnplacedMultiUnion const *>(solid)) {
    for (size_t i = 0; i < multi->GetNumberOfSolids(); ++i) {
      auto const *placed = multi->GetNode(i);
      Transformation3D transform(frame);
      transform.MultiplyFromRight(*placed->GetTransformation());
      CollectUnion(placed->GetUnplacedVolume(), transform, components);
    }
  } else {
    components.push_back({solid, frame});
  }
}

bool CanBuildUnion(std::vector<Component> const &components, size_t minimum)
{
  if (components.size() < minimum || components.size() < 2) return false;
  for (auto const &component : components)
    if (!BooleanHelper::IsFiniteBody(component.solid)) return false;
  return true;
}

VUnplacedVolume const *MakeBoolean(BooleanOperation op, VPlacedVolume const *left, VPlacedVolume const *right)
{
  switch (op) {
  case kUnion:
    return new UnplacedBooleanVolume<kUnion>(op, left, right);
  case kIntersection:
    return new UnplacedBooleanVolume<kIntersection>(op, left, right);
  case kSubtraction:
    return new UnplacedBooleanVolume<kSubtraction>(op, left, right);
  }
  return nullptr;
}

} // namespace

VPlacedVolume const *BooleanFactory::Place(VUnplacedVolume const *solid, Transformation3D const &transform)
{
  auto &logical = fLogicalVolumes[solid];
  if (!logical) logical = new LogicalVolume(solid);
  return logical->Place(new Transformation3D(transform));
}

VPlacedVolume const *BooleanFactory::ConvertPlaced(VPlacedVolume const *placed)
{
  auto const *solid = Create(placed->GetUnplacedVolume());
  return solid == placed->GetUnplacedVolume() ? placed : Place(solid, *placed->GetTransformation());
}

VUnplacedVolume const *BooleanFactory::Create(VUnplacedVolume const *solid)
{
  if (!solid || fMinComponents == 0) return solid;
  auto found = fResults.find(solid);
  if (found != fResults.end()) return found->second;

  VUnplacedVolume const *result = solid;
  auto const *boolean           = BooleanHelper::GetBooleanStruct(solid);
  if (boolean) {
    std::vector<Component> components;
    VUnplacedVolume const *body = solid;
    Transformation3D bodyFrame;
    auto const *leftBoolean    = BooleanHelper::GetBooleanStruct(boolean->fLeftVolume->GetUnplacedVolume());
    const bool alreadyCombined = boolean->fOp == kSubtraction &&
                                 dynamic_cast<UnplacedMultiUnion const *>(boolean->fRightVolume->GetUnplacedVolume()) &&
                                 (!leftBoolean || leftBoolean->fOp != kSubtraction);
    if (boolean->fOp == kUnion) {
      CollectUnion(solid, bodyFrame, components);
    } else if (boolean->fOp == kSubtraction && !alreadyCombined) {
      // (A-B)-C = A-(B union C). Accumulate every placement on the left
      // spine before expressing its right cutter in the outer Boolean frame.
      auto const *node = boolean;
      do {
        Transformation3D cutterFrame(bodyFrame);
        cutterFrame.MultiplyFromRight(*node->fRightVolume->GetTransformation());
        CollectUnion(node->fRightVolume->GetUnplacedVolume(), cutterFrame, components);
        bodyFrame.MultiplyFromRight(*node->fLeftVolume->GetTransformation());
        body = node->fLeftVolume->GetUnplacedVolume();
        node = BooleanHelper::GetBooleanStruct(body);
      } while (node && node->fOp == kSubtraction);
    }

    if (CanBuildUnion(components, fMinComponents)) {
      auto *multi = new UnplacedMultiUnion();
      for (auto const &component : components)
        multi->AddNode(Place(Create(component.solid), component.transform));
      multi->Close();
      result = boolean->fOp == kUnion ? static_cast<VUnplacedVolume const *>(multi)
                                      : MakeBoolean(kSubtraction, Place(Create(body), bodyFrame),
                                                    Place(multi, Transformation3D::kIdentity));
    } else {
      auto const *left  = ConvertPlaced(boolean->fLeftVolume);
      auto const *right = ConvertPlaced(boolean->fRightVolume);
      if (left != boolean->fLeftVolume || right != boolean->fRightVolume)
        result = MakeBoolean(boolean->fOp, left, right);
    }
  } else if (auto const *scaled = dynamic_cast<UnplacedScaledShape const *>(solid)) {
    auto const *placed = ConvertPlaced(scaled->GetPlaced());
    if (placed != scaled->GetPlaced()) {
      auto const &scale = scaled->GetScale().Scale();
      result            = new UnplacedScaledShape(placed, scale.x(), scale.y(), scale.z());
    }
  } else if (auto const *multi = dynamic_cast<UnplacedMultiUnion const *>(solid)) {
    std::vector<VPlacedVolume const *> components;
    bool changed = false;
    for (size_t i = 0; i < multi->GetNumberOfSolids(); ++i) {
      auto const *placed = ConvertPlaced(multi->GetNode(i));
      components.push_back(placed);
      changed |= placed != multi->GetNode(i);
    }
    if (changed) {
      auto *replacement = new UnplacedMultiUnion();
      for (auto const *placed : components)
        replacement->AddNode(placed);
      replacement->Close();
      result = replacement;
    }
  }

  fResults[solid]  = result;
  fResults[result] = result;
  return result;
}

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

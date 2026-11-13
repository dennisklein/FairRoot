//===-- source/BVHSubdivisionRecipe.cpp ---------------------------*- C++ -*-===//
//
/// @file BVHSubdivisionRecipe.cpp
/// @brief Strict binary reader and writer for persisted BVH split decisions.
//
//===----------------------------------------------------------------------===//

#include "BVHSubdivisionRecipe.h"
#include "VecGeom/management/BVHSubdivisionOptions.h"

#include <array>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <unordered_set>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
namespace bvhrecipe {
namespace {

// The checksum covers every byte before its final U64 value. All integer fields
// use the BinaryReader/BinaryWriter little-endian helpers; no raw structs or
// floating-point coordinates are part of the recipe file format.
constexpr std::array<unsigned char, 16> kMagic = {'V', 'E', 'C', 'G', 'E', 'O', 'M', '_',
                                                  'B', 'V', 'H', '_', 'S', 'P', 'L', 'T'};
constexpr std::uint64_t kChecksumOffset        = UINT64_C(14695981039346656037);
constexpr std::uint64_t kChecksumPrime         = UINT64_C(1099511628211);

void UpdateChecksum(std::uint64_t &checksum, void const *data, std::size_t size)
{
  auto const *bytes = static_cast<unsigned char const *>(data);
  for (std::size_t index = 0; index < size; ++index) {
    checksum ^= bytes[index];
    checksum *= kChecksumPrime;
  }
}

bool Fail(std::string message, std::string &error)
{
  error = std::move(message);
  return false;
}

class BinaryWriter {
private:
  std::ofstream fOutput;
  std::uint64_t fChecksum{kChecksumOffset};

  void RawBytes(void const *data, std::size_t size) { fOutput.write(static_cast<char const *>(data), size); }

public:
  explicit BinaryWriter(std::string const &fileName) : fOutput(fileName, std::ios::binary | std::ios::trunc) {}

  bool IsOpen() const { return fOutput.is_open(); }
  bool Finish()
  {
    unsigned char checksumBytes[8];
    auto checksum = fChecksum;
    for (int index = 0; index < 8; ++index) {
      checksumBytes[index] = static_cast<unsigned char>(checksum & UINT64_C(0xff));
      checksum >>= 8;
    }
    RawBytes(checksumBytes, sizeof(checksumBytes));
    fOutput.flush();
    if (!fOutput) return false;
    fOutput.close();
    return !fOutput.fail();
  }

  void Bytes(void const *data, std::size_t size)
  {
    RawBytes(data, size);
    UpdateChecksum(fChecksum, data, size);
  }

  void U32(std::uint32_t value)
  {
    unsigned char bytes[4];
    for (int index = 0; index < 4; ++index) {
      bytes[index] = static_cast<unsigned char>(value & UINT32_C(0xff));
      value >>= 8;
    }
    Bytes(bytes, sizeof(bytes));
  }

  void U8(std::uint8_t value) { Bytes(&value, sizeof(value)); }

  void U64(std::uint64_t value)
  {
    unsigned char bytes[8];
    for (int index = 0; index < 8; ++index) {
      bytes[index] = static_cast<unsigned char>(value & UINT64_C(0xff));
      value >>= 8;
    }
    Bytes(bytes, sizeof(bytes));
  }

  void String(std::string const &value)
  {
    U64(value.size());
    Bytes(value.data(), value.size());
  }
};

class BinaryReader {
private:
  std::ifstream fInput;
  std::size_t fRemaining{0};
  std::uint64_t fChecksum{kChecksumOffset};

  bool RawBytes(void *data, std::size_t size)
  {
    if (size > fRemaining) return false;
    fInput.read(static_cast<char *>(data), size);
    if (!fInput) return false;
    fRemaining -= size;
    return true;
  }

public:
  BinaryReader(std::string const &fileName, std::size_t size) : fInput(fileName, std::ios::binary), fRemaining(size) {}

  bool IsOpen() const { return fInput.is_open(); }
  std::size_t Remaining() const { return fRemaining; }

  bool Bytes(void *data, std::size_t size)
  {
    if (!RawBytes(data, size)) return false;
    UpdateChecksum(fChecksum, data, size);
    return true;
  }

  bool ChecksumMatches()
  {
    if (fRemaining != sizeof(std::uint64_t)) return false;
    unsigned char bytes[8];
    if (!RawBytes(bytes, sizeof(bytes))) return false;
    std::uint64_t stored = 0;
    for (int index = 7; index >= 0; --index)
      stored = (stored << 8) | bytes[index];
    return stored == fChecksum;
  }

  bool U32(std::uint32_t &value)
  {
    unsigned char bytes[4];
    if (!Bytes(bytes, sizeof(bytes))) return false;
    value = 0;
    for (int index = 3; index >= 0; --index)
      value = (value << 8) | bytes[index];
    return true;
  }

  bool U8(std::uint8_t &value) { return Bytes(&value, sizeof(value)); }

  bool U64(std::uint64_t &value)
  {
    unsigned char bytes[8];
    if (!Bytes(bytes, sizeof(bytes))) return false;
    value = 0;
    for (int index = 7; index >= 0; --index)
      value = (value << 8) | bytes[index];
    return true;
  }

  bool String(std::string &value)
  {
    std::uint64_t size = 0;
    if (!U64(size) || size > fRemaining || size > std::numeric_limits<std::size_t>::max()) return false;
    value.resize(static_cast<std::size_t>(size));
    return Bytes(value.data(), value.size());
  }
};

bool FileSize(std::string const &fileName, std::size_t &size)
{
  std::ifstream input(fileName, std::ios::binary | std::ios::ate);
  if (!input) return false;
  auto const end = input.tellg();
  if (end < 0 || static_cast<unsigned long long>(end) > std::numeric_limits<std::size_t>::max()) return false;
  size = static_cast<std::size_t>(end);
  return true;
}

bool ReadCount(BinaryReader &input, std::size_t &count, std::size_t minimumElementBytes = 1)
{
  std::uint64_t value = 0;
  if (!input.U64(value) || value > std::numeric_limits<std::size_t>::max()) return false;
  if (minimumElementBytes > 0 && value > input.Remaining() / minimumElementBytes) return false;
  count = static_cast<std::size_t>(value);
  return true;
}

bool ValidateForWrite(Recipe const &recipe, std::string &error)
{
  if (recipe.precisionBytes <= 0 || recipe.bvhRealBytes <= 0)
    return Fail("cannot write invalid recipe precision metadata", error);

  std::unordered_set<int> logicalIds;
  for (auto const &volume : recipe.volumes) {
    if (volume.logicalId < 0 || !logicalIds.insert(volume.logicalId).second ||
        volume.daughterCount > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
      return Fail("cannot write invalid or duplicate logical-volume entry", error);

    int previousDaughterId = -1;
    for (auto const &daughter : volume.splitDaughters) {
      if (daughter.daughterId <= previousDaughterId ||
          static_cast<std::size_t>(daughter.daughterId) >= volume.daughterCount || daughter.axis < 0 ||
          daughter.axis > 2 || daughter.M <= 1 || daughter.M > kMaximumAABBSubdivisionM)
        return Fail("cannot write invalid split decision", error);
      previousDaughterId = daughter.daughterId;
    }
  }
  return true;
}

std::filesystem::path TemporaryPublishPath(std::filesystem::path const &target)
{
  static std::atomic<std::uint64_t> sequence{0};
  auto const timestamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::uint64_t nonce  = reinterpret_cast<std::uintptr_t>(&sequence);
  try {
    std::random_device entropy;
    nonce ^= (static_cast<std::uint64_t>(entropy()) << 32) ^ static_cast<std::uint64_t>(entropy());
  } catch (...) {
    // Timestamp, process address, and sequence still make a collision remote
    // when a standard-library entropy source is unavailable.
  }
  auto path = target;
  path +=
      ".tmp." + std::to_string(timestamp) + "." + std::to_string(nonce) + "." + std::to_string(sequence.fetch_add(1));
  return path;
}

bool DestinationExists(std::filesystem::path const &path)
{
  std::error_code ignored;
  return std::filesystem::exists(path, ignored);
}

} // namespace

bool Read(std::string const &fileName, Recipe &recipe, std::string &error)
{
  std::size_t fileBytes = 0;
  if (!FileSize(fileName, fileBytes)) return Fail("cannot open recipe file '" + fileName + "'", error);
  if (fileBytes < kMagic.size() + sizeof(std::uint32_t) + sizeof(std::uint64_t))
    return Fail("recipe file is truncated", error);

  BinaryReader input(fileName, fileBytes);
  std::array<unsigned char, kMagic.size()> magic;
  std::uint32_t formatVersion = 0;
  if (!input.IsOpen() || !input.Bytes(magic.data(), magic.size()) || magic != kMagic || !input.U32(formatVersion))
    return Fail("not a VecGeom BVH AABB subdivision recipe", error);
  if (formatVersion != kFormatVersion)
    return Fail("unsupported recipe format version " + std::to_string(formatVersion), error);

  Recipe parsed;
  std::uint32_t integer = 0;
  if (!input.String(parsed.vecgeomVersion) || !input.U32(integer) || integer > std::numeric_limits<int>::max())
    return Fail("invalid recipe metadata", error);
  parsed.precisionBytes = static_cast<int>(integer);
  if (!input.U32(integer) || integer > std::numeric_limits<int>::max()) return Fail("invalid recipe metadata", error);
  parsed.bvhRealBytes = static_cast<int>(integer);
  if (!input.String(parsed.worldLabel) || !input.U64(parsed.geometryFingerprint))
    return Fail("invalid recipe metadata", error);

  std::size_t volumeCount = 0;
  if (!ReadCount(input, volumeCount, 36)) return Fail("invalid volume count", error);
  parsed.volumes.reserve(volumeCount);
  std::unordered_set<int> logicalIds;
  for (std::size_t volumeIndex = 0; volumeIndex < volumeCount; ++volumeIndex) {
    Volume volume;
    if (!input.U32(integer) || integer > std::numeric_limits<int>::max() ||
        !logicalIds.insert(static_cast<int>(integer)).second || !input.String(volume.label) ||
        !input.U64(volume.fingerprint))
      return Fail("invalid or duplicate logical-volume entry", error);
    volume.logicalId = static_cast<int>(integer);

    if (!ReadCount(input, volume.daughterCount, 0) ||
        volume.daughterCount > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
      return Fail("invalid daughter count", error);
    std::size_t splitCount = 0;
    if (!ReadCount(input, splitCount, sizeof(std::uint32_t) + sizeof(std::uint8_t) + sizeof(std::uint32_t)) ||
        splitCount > volume.daughterCount)
      return Fail("invalid split-daughter count", error);
    volume.splitDaughters.reserve(splitCount);
    int previousDaughterId = -1;
    for (std::size_t splitIndex = 0; splitIndex < splitCount; ++splitIndex) {
      Daughter daughter;
      std::uint8_t axis = 0;
      std::uint32_t M   = 0;
      if (!input.U32(integer) || integer > static_cast<std::uint32_t>(std::numeric_limits<int>::max()) ||
          static_cast<int>(integer) <= previousDaughterId || integer >= volume.daughterCount || !input.U8(axis) ||
          axis > 2 || !input.U32(M) || M <= 1 || M > static_cast<std::uint32_t>(kMaximumAABBSubdivisionM))
        return Fail("invalid or unordered split decision", error);
      previousDaughterId  = static_cast<int>(integer);
      daughter.daughterId = static_cast<int>(integer);
      daughter.axis       = axis;
      daughter.M          = static_cast<int>(M);
      volume.splitDaughters.push_back(std::move(daughter));
    }
    parsed.volumes.push_back(std::move(volume));
  }
  if (!input.ChecksumMatches()) return Fail("recipe payload checksum mismatch or trailing data", error);

  recipe = std::move(parsed);
  error.clear();
  return true;
}

bool Write(std::string const &fileName, Recipe const &recipe, std::string &error)
{
  if (!ValidateForWrite(recipe, error)) return false;

  BinaryWriter output(fileName);
  if (!output.IsOpen()) return Fail("cannot create recipe file '" + fileName + "'", error);

  output.Bytes(kMagic.data(), kMagic.size());
  output.U32(kFormatVersion);
  output.String(recipe.vecgeomVersion);
  output.U32(static_cast<std::uint32_t>(recipe.precisionBytes));
  output.U32(static_cast<std::uint32_t>(recipe.bvhRealBytes));
  output.String(recipe.worldLabel);
  output.U64(recipe.geometryFingerprint);
  output.U64(recipe.volumes.size());
  for (auto const &volume : recipe.volumes) {
    output.U32(static_cast<std::uint32_t>(volume.logicalId));
    output.String(volume.label);
    output.U64(volume.fingerprint);
    output.U64(volume.daughterCount);
    output.U64(volume.splitDaughters.size());
    for (auto const &daughter : volume.splitDaughters) {
      output.U32(static_cast<std::uint32_t>(daughter.daughterId));
      output.U8(static_cast<std::uint8_t>(daughter.axis));
      output.U32(static_cast<std::uint32_t>(daughter.M));
    }
  }
  if (!output.Finish()) return Fail("failed while writing recipe file '" + fileName + "'", error);

  error.clear();
  return true;
}

PublishResult WriteIfAbsent(std::string const &fileName, Recipe const &recipe, std::string &error)
{
  std::filesystem::path const target(fileName);
  std::error_code statusError;
  if (std::filesystem::exists(target, statusError)) {
    error.clear();
    return PublishResult::kExists;
  }
  if (statusError) {
    error = "cannot inspect recipe destination '" + fileName + "': " + statusError.message();
    return PublishResult::kError;
  }

  auto const temporary = TemporaryPublishPath(target);
  if (!Write(temporary.string(), recipe, error)) {
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    return PublishResult::kError;
  }

  // A hard link publishes the fully written payload atomically and fails when
  // another process has already won. Do not fall back to copying: a concurrent
  // reader could observe a partially copied destination.
  std::error_code publishError;
  std::filesystem::create_hard_link(temporary, target, publishError);

  std::error_code removeError;
  std::filesystem::remove(temporary, removeError);
  if (!publishError) {
    error.clear();
    return PublishResult::kCreated;
  }
  if (DestinationExists(target)) {
    error.clear();
    return PublishResult::kExists;
  }

  error = "cannot publish recipe file '" + fileName + "': " + publishError.message();
  return PublishResult::kError;
}

} // namespace bvhrecipe
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

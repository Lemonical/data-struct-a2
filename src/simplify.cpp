#include "simplify.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

#include "geometry.hpp"
#include "validation.hpp"

namespace {

// This epsilon keeps tie-break comparisons stable when doubles are nearly equal
constexpr double kEpsilon = 1e-12;

// This tolerance controls ring-area restoration stop criteria in relative terms
constexpr double kAreaRestoreRelativeTolerance = 1e-9;

// This tolerance guards very small absolute area differences near zero-area limits
constexpr double kAreaRestoreAbsoluteTolerance = 1e-9;

// This cap keeps the restoration pass deterministic and bounded
constexpr std::size_t kMaxAreaRestoreIterations = 12U;

// This captures one removable-vertex option with deterministic ordering fields
struct RemovalCandidate {
    std::size_t ringIndex = 0U;
    std::size_t vertexIndex = 0U;
    int ringId = -1;
    atpps::Point point;
    double shapeError = std::numeric_limits<double>::infinity();
    double areaDelta = std::numeric_limits<double>::infinity();
};

// This computes squared distance from p to segment a-b for local shape-error scoring
double SquaredDistanceToSegment(const atpps::Point& p, const atpps::Point& a, const atpps::Point& b) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double lengthSquared = (dx * dx) + (dy * dy);

    if (lengthSquared <= kEpsilon) {
        const double px = p.x - a.x;
        const double py = p.y - a.y;
        return (px * px) + (py * py);
    }

    const double projection = ((p.x - a.x) * dx + (p.y - a.y) * dy) / lengthSquared;
    const double t = (projection < 0.0) ? 0.0 : ((projection > 1.0) ? 1.0 : projection);
    const double nearestX = a.x + (t * dx);
    const double nearestY = a.y + (t * dy);
    const double diffX = p.x - nearestX;
    const double diffY = p.y - nearestY;
    return (diffX * diffX) + (diffY * diffY);
}

// This returns the signed area of triangle a-b-c used for local area-change scoring
double SignedTriangleArea(const atpps::Point& a, const atpps::Point& b, const atpps::Point& c) {
    const double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    return 0.5 * cross;
}

// This builds deterministic local metrics for one candidate removal
RemovalCandidate BuildCandidate(
    const atpps::Ring& ring,
    const std::size_t ringIndex,
    const std::size_t vertexIndex
) {
    const std::size_t count = ring.vertices.size();
    const std::size_t prevIndex = (vertexIndex + count - 1U) % count;
    const std::size_t nextIndex = (vertexIndex + 1U) % count;

    const atpps::Point& prev = ring.vertices[prevIndex];
    const atpps::Point& current = ring.vertices[vertexIndex];
    const atpps::Point& next = ring.vertices[nextIndex];

    RemovalCandidate candidate;
    candidate.ringIndex = ringIndex;
    candidate.vertexIndex = vertexIndex;
    candidate.ringId = ring.ringId;
    candidate.point = current;

    // This favors removals that keep the local curve closest to its neighbor chord
    candidate.shapeError = SquaredDistanceToSegment(current, prev, next);

    // This favors removals that perturb local signed area as little as possible
    candidate.areaDelta = std::abs(SignedTriangleArea(prev, current, next));
    return candidate;
}

// This compares doubles with epsilon while keeping strict weak ordering behavior
bool LessWithTolerance(const double lhs, const double rhs) {
    if (lhs + kEpsilon < rhs) {
        return true;
    }
    return false;
}

// This enforces deterministic candidate ordering for global best-choice selection
bool IsBetterCandidate(const RemovalCandidate& lhs, const RemovalCandidate& rhs) {
    if (LessWithTolerance(lhs.shapeError, rhs.shapeError)) {
        return true;
    }
    if (LessWithTolerance(rhs.shapeError, lhs.shapeError)) {
        return false;
    }

    if (LessWithTolerance(lhs.areaDelta, rhs.areaDelta)) {
        return true;
    }
    if (LessWithTolerance(rhs.areaDelta, lhs.areaDelta)) {
        return false;
    }

    if (lhs.ringId != rhs.ringId) {
        return lhs.ringId < rhs.ringId;
    }

    if (lhs.vertexIndex != rhs.vertexIndex) {
        return lhs.vertexIndex < rhs.vertexIndex;
    }

    if (LessWithTolerance(lhs.point.x, rhs.point.x)) {
        return true;
    }
    if (LessWithTolerance(rhs.point.x, lhs.point.x)) {
        return false;
    }

    if (LessWithTolerance(lhs.point.y, rhs.point.y)) {
        return true;
    }
    if (LessWithTolerance(rhs.point.y, lhs.point.y)) {
        return false;
    }

    return false;
}

// This applies one vertex deletion by index in one ring in a copied polygon state
atpps::Polygon ApplyRemoval(const atpps::Polygon& polygon, const RemovalCandidate& candidate) {
    atpps::Polygon edited = polygon;
    auto& ring = edited.rings[candidate.ringIndex];
    ring.vertices.erase(ring.vertices.begin() + static_cast<std::ptrdiff_t>(candidate.vertexIndex));
    return edited;
}

// This scans all legal removals and returns the best one that preserves topology
bool FindBestLegalRemoval(const atpps::Polygon& polygon, RemovalCandidate& bestCandidate) {
    bool foundAny = false;

    for (std::size_t ringIndex = 0U; ringIndex < polygon.rings.size(); ++ringIndex) {
        const atpps::Ring& ring = polygon.rings[ringIndex];

        // This keeps rings valid by never dropping below a closed triangle
        if (ring.vertices.size() <= 3U) {
            continue;
        }

        for (std::size_t vertexIndex = 0U; vertexIndex < ring.vertices.size(); ++vertexIndex) {
            const RemovalCandidate candidate = BuildCandidate(ring, ringIndex, vertexIndex);
            const atpps::Polygon edited = ApplyRemoval(polygon, candidate);

            std::string topologyError;
            if (!atpps::ValidatePolygonTopology(edited, topologyError)) {
                continue;
            }

            if (!foundAny || IsBetterCandidate(candidate, bestCandidate)) {
                bestCandidate = candidate;
                foundAny = true;
            }
        }
    }

    return foundAny;
}

// This computes arithmetic-mean center used as a stable scaling pivot for one ring
atpps::Point ComputeRingCentroid(const atpps::Ring& ring) {
    atpps::Point centroid;
    if (ring.vertices.empty()) {
        return centroid;
    }

    for (const atpps::Point& point : ring.vertices) {
        centroid.x += point.x;
        centroid.y += point.y;
    }

    const double invCount = 1.0 / static_cast<double>(ring.vertices.size());
    centroid.x *= invCount;
    centroid.y *= invCount;
    return centroid;
}

// This returns a ring scaled around a center so signed area changes while local shape stays similar
atpps::Ring ScaleRingAroundCenter(const atpps::Ring& ring, const atpps::Point& center, const double scale) {
    atpps::Ring scaled = ring;
    for (atpps::Point& point : scaled.vertices) {
        point.x = center.x + (point.x - center.x) * scale;
        point.y = center.y + (point.y - center.y) * scale;
    }
    return scaled;
}

// This checks whether two signed-area values are close enough under mixed relative and absolute tolerance
bool AreasCloseEnough(const double lhs, const double rhs) {
    const double delta = std::abs(lhs - rhs);
    const double scale = std::max(std::abs(lhs), std::max(std::abs(rhs), 1.0));
    return delta <= (kAreaRestoreAbsoluteTolerance + kAreaRestoreRelativeTolerance * scale);
}

// This appends one note fragment while preserving a single-line stderr note style
void AppendNoteMessage(std::string& noteMessage, const std::string& fragment) {
    if (fragment.empty()) {
        return;
    }

    if (!noteMessage.empty()) {
        noteMessage += " | ";
    }
    noteMessage += fragment;
}

// This attempts deterministic safe ring-area restoration by iterative damped scaling with topology checks
bool RestoreRingAreaSafely(
    atpps::Polygon& polygon,
    const std::size_t ringIndex,
    const double targetSignedArea,
    std::string& failureReason
) {
    atpps::Ring& ring = polygon.rings[ringIndex];

    // This keeps restoration logic safe for non-polygonal rings
    if (ring.vertices.size() < 3U) {
        failureReason = "ring has fewer than 3 vertices";
        return false;
    }

    // This uses one stable pivot so repeated scaling stays deterministic
    const atpps::Point center = ComputeRingCentroid(ring);

    // This uses deterministic damping steps from full correction toward conservative updates
    const double dampingFactors[] = {1.0, 0.5, 0.25, 0.125, 0.0625};

    for (std::size_t iteration = 0U; iteration < kMaxAreaRestoreIterations; ++iteration) {
        const double currentArea = atpps::ComputeSignedArea(ring);
        if (AreasCloseEnough(currentArea, targetSignedArea)) {
            failureReason.clear();
            return true;
        }

        // This avoids unstable correction when area sign is incompatible with target sign
        if (currentArea * targetSignedArea <= 0.0) {
            failureReason = "current and target ring signed areas have incompatible signs";
            return false;
        }

        // This avoids division by very small area values that would explode scaling factors
        if (std::abs(currentArea) <= kEpsilon) {
            failureReason = "current ring area is too close to zero for stable scaling";
            return false;
        }

        const double idealScale = std::sqrt(std::abs(targetSignedArea / currentArea));

        bool foundImprovement = false;
        atpps::Ring bestRing = ring;
        double bestError = std::abs(currentArea - targetSignedArea);

        for (const double damping : dampingFactors) {
            const double candidateScale = 1.0 + damping * (idealScale - 1.0);

            // This rejects collapsed or flipped scaling factors before geometry checks
            if (candidateScale <= kEpsilon) {
                continue;
            }

            atpps::Polygon candidatePolygon = polygon;
            candidatePolygon.rings[ringIndex] = ScaleRingAroundCenter(ring, center, candidateScale);

            std::string topologyError;
            if (!atpps::ValidatePolygonTopology(candidatePolygon, topologyError)) {
                continue;
            }

            const double candidateArea = atpps::ComputeSignedArea(candidatePolygon.rings[ringIndex]);
            const double candidateError = std::abs(candidateArea - targetSignedArea);

            // This accepts only strict improvement so iteration progress is monotonic
            if (candidateError + kEpsilon < bestError) {
                bestError = candidateError;
                bestRing = candidatePolygon.rings[ringIndex];
                foundImprovement = true;
            }
        }

        // This reports constrained restoration when no safe improving step exists
        if (!foundImprovement) {
            failureReason = "no topology-safe scaling improvement found";
            return false;
        }

        ring = bestRing;
    }

    const double finalArea = atpps::ComputeSignedArea(ring);
    if (AreasCloseEnough(finalArea, targetSignedArea)) {
        failureReason.clear();
        return true;
    }

    failureReason = "iteration cap reached before meeting area tolerance";
    return false;
}

// This restores ring signed areas toward input targets and records constrained-stop diagnostics
void RestorePolygonRingAreas(const atpps::Polygon& inputPolygon, atpps::Polygon& simplifiedPolygon, std::string& noteMessage) {
    std::map<int, double> targetAreaByRingId;
    for (const atpps::Ring& ring : inputPolygon.rings) {
        targetAreaByRingId[ring.ringId] = atpps::ComputeSignedArea(ring);
    }

    for (std::size_t ringIndex = 0U; ringIndex < simplifiedPolygon.rings.size(); ++ringIndex) {
        const int ringId = simplifiedPolygon.rings[ringIndex].ringId;
        const auto targetAreaIt = targetAreaByRingId.find(ringId);
        if (targetAreaIt == targetAreaByRingId.end()) {
            AppendNoteMessage(noteMessage, "ring " + std::to_string(ringId) + " has no input area target");
            continue;
        }

        const double targetArea = targetAreaIt->second;
        const double currentArea = atpps::ComputeSignedArea(simplifiedPolygon.rings[ringIndex]);
        if (AreasCloseEnough(currentArea, targetArea)) {
            continue;
        }

        std::string failureReason;
        const bool restored = RestoreRingAreaSafely(simplifiedPolygon, ringIndex, targetArea, failureReason);
        if (!restored) {
            std::ostringstream message;
            message << "ring " << ringId << " area restoration constrained: " << failureReason;
            AppendNoteMessage(noteMessage, message.str());
        }
    }
}

}  // namespace

namespace atpps {

SimplificationResult SimplifyPolygonToTarget(
    const Polygon& inputPolygon,
    const std::size_t targetVertices,
    std::string& noteMessage
) {
    // This starts from input and performs deterministic legal removals toward target
    SimplificationResult result;
    result.polygon = inputPolygon;
    result.requestedTarget = targetVertices;
    noteMessage.clear();

    std::size_t currentVertices = CountTotalVertices(result.polygon);
    if (targetVertices >= currentVertices) {
        result.finalVertexCount = currentVertices;
        result.reachedExactTarget = (result.finalVertexCount == targetVertices);
        if (!result.reachedExactTarget) {
            noteMessage = "target is greater than current vertex count so no removals were needed";
        }
        return result;
    }

    while (currentVertices > targetVertices) {
        RemovalCandidate bestCandidate;
        const bool foundRemoval = FindBestLegalRemoval(result.polygon, bestCandidate);
        if (!foundRemoval) {
            break;
        }

        result.polygon = ApplyRemoval(result.polygon, bestCandidate);
        currentVertices = CountTotalVertices(result.polygon);
    }

    // This performs post-adjustment to restore per-ring signed areas when safe
    RestorePolygonRingAreas(inputPolygon, result.polygon, noteMessage);

    result.finalVertexCount = currentVertices;
    result.reachedExactTarget = (result.finalVertexCount == targetVertices);

    // This reports constrained-stop cases through stderr while preserving stdout format
    if (!result.reachedExactTarget) {
        noteMessage = "target could not be reached under topology constraints; stopped at feasible vertex count";
    }

    return result;
}

}  // namespace atpps

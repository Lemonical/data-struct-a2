#include "simplify.hpp"

#include <cmath>
#include <limits>

#include "geometry.hpp"
#include "validation.hpp"

namespace {

// This epsilon keeps tie-break comparisons stable when doubles are nearly equal
constexpr double kEpsilon = 1e-12;

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

    result.finalVertexCount = currentVertices;
    result.reachedExactTarget = (result.finalVertexCount == targetVertices);

    // This reports constrained-stop cases through stderr while preserving stdout format
    if (!result.reachedExactTarget) {
        noteMessage = "target could not be reached under topology constraints; stopped at feasible vertex count";
    }

    return result;
}

}  // namespace atpps

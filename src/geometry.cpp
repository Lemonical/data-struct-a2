#include "geometry.hpp"

#include <cmath>
#include <map>

namespace atpps {

double ComputeSignedArea(const Ring& ring) {
    // This guards against degenerate rings so stats stay predictable
    if (ring.vertices.size() < 3U) {
        return 0.0;
    }

    // This uses the shoelace sum with implicit edge from last back to first
    double twiceArea = 0.0;
    const std::size_t count = ring.vertices.size();
    for (std::size_t i = 0; i < count; ++i) {
        const Point& current = ring.vertices[i];
        const Point& next = ring.vertices[(i + 1U) % count];
        twiceArea += (current.x * next.y) - (next.x * current.y);
    }

    // This converts from twice-area to true area while preserving sign
    return 0.5 * twiceArea;
}

double ComputeTotalSignedArea(const Polygon& polygon) {
    // This folds per-ring signed area into one total metric
    double total = 0.0;
    for (const Ring& ring : polygon.rings) {
        total += ComputeSignedArea(ring);
    }
    return total;
}

double ComputeTotalArealDisplacement(const Polygon& inputPolygon, const Polygon& outputPolygon) {
    // This maps ring id to signed area so we can compare matching rings deterministically
    std::map<int, double> inputAreaByRing;
    for (const Ring& ring : inputPolygon.rings) {
        inputAreaByRing[ring.ringId] = ComputeSignedArea(ring);
    }

    // This maps output ring areas with the same ring id key space
    std::map<int, double> outputAreaByRing;
    for (const Ring& ring : outputPolygon.rings) {
        outputAreaByRing[ring.ringId] = ComputeSignedArea(ring);
    }

    // This accumulates absolute area delta per ring as a stable baseline metric
    double totalDisplacement = 0.0;
    for (const auto& [ringId, inputArea] : inputAreaByRing) {
        const auto outputIt = outputAreaByRing.find(ringId);
        const double outputArea = (outputIt != outputAreaByRing.end()) ? outputIt->second : 0.0;
        totalDisplacement += std::abs(inputArea - outputArea);
    }

    // This also accounts for rings that only exist in output data
    for (const auto& [ringId, outputArea] : outputAreaByRing) {
        if (inputAreaByRing.find(ringId) == inputAreaByRing.end()) {
            totalDisplacement += std::abs(outputArea);
        }
    }

    return totalDisplacement;
}

std::size_t CountTotalVertices(const Polygon& polygon) {
    // This sums all ring vertex counts for global target handling
    std::size_t total = 0U;
    for (const Ring& ring : polygon.rings) {
        total += ring.vertices.size();
    }
    return total;
}

}  // namespace atpps

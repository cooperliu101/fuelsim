#include "spatial_common.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace fuelsim::spatial_detail {
namespace {
bool farther_contact_search_entry(const ContactSearchQueueEntry& first, const ContactSearchQueueEntry& second) {
    return first.distance_squared > second.distance_squared;
}
} // namespace

void ContactSearchTree::build(std::vector<ContactSearchBox> boxes) {
    _boxes = std::move(boxes);
    _nodes.clear();
    if (_boxes.empty()) return;
    std::vector<std::size_t> indices(_boxes.size());
    std::iota(indices.begin(), indices.end(), 0);
    _nodes.reserve(2 * _boxes.size() - 1);
    (void)build_node(indices, 0, indices.size());
}

void ContactSearchTree::refit(const std::vector<ContactSearchBox>& boxes) {
    if (!can_refit(boxes.size())) throw std::invalid_argument("Contact search tree refit layout mismatch");
    for (std::size_t offset = _nodes.size(); offset-- > 0;) {
        Node& node = _nodes[offset];
        if (node.item != std::numeric_limits<std::size_t>::max()) {
            if (node.item >= boxes.size() || boxes[node.item].item != node.item)
                throw std::invalid_argument("Contact search tree refit items must use dense indices");
            node.minimum = boxes[node.item].minimum;
            node.maximum = boxes[node.item].maximum;
            continue;
        }
        for (std::size_t component = 0; component < 3; ++component) {
            node.minimum[component] =
                std::min(_nodes[node.left].minimum[component], _nodes[node.right].minimum[component]);
            node.maximum[component] =
                std::max(_nodes[node.left].maximum[component], _nodes[node.right].maximum[component]);
        }
    }
}

std::size_t ContactSearchTree::build_node(std::vector<std::size_t>& indices, std::size_t begin, std::size_t end) {
    const std::size_t node_index = _nodes.size();
    _nodes.push_back({});
    Node bounds;
    bounds.minimum.fill(std::numeric_limits<double>::infinity());
    bounds.maximum.fill(-std::numeric_limits<double>::infinity());
    for (std::size_t offset = begin; offset < end; ++offset)
        for (std::size_t component = 0; component < 3; ++component) {
            bounds.minimum[component] = std::min(bounds.minimum[component], _boxes[indices[offset]].minimum[component]);
            bounds.maximum[component] = std::max(bounds.maximum[component], _boxes[indices[offset]].maximum[component]);
        }
    if (end - begin == 1) {
        bounds.item = _boxes[indices[begin]].item;
        _nodes[node_index] = bounds;
        return node_index;
    }
    std::size_t axis = 0;
    for (std::size_t component = 1; component < 3; ++component)
        if (bounds.maximum[component] - bounds.minimum[component] > bounds.maximum[axis] - bounds.minimum[axis])
            axis = component;
    const std::size_t middle = begin + (end - begin) / 2;
    std::nth_element(indices.begin() + static_cast<std::ptrdiff_t>(begin),
        indices.begin() + static_cast<std::ptrdiff_t>(middle), indices.begin() + static_cast<std::ptrdiff_t>(end),
        [this, axis](std::size_t first, std::size_t second) {
            const double first_center = 0.5 * (_boxes[first].minimum[axis] + _boxes[first].maximum[axis]);
            const double second_center = 0.5 * (_boxes[second].minimum[axis] + _boxes[second].maximum[axis]);
            return first_center < second_center ||
                   (first_center == second_center && _boxes[first].item < _boxes[second].item);
        });
    bounds.left = build_node(indices, begin, middle);
    bounds.right = build_node(indices, middle, end);
    _nodes[node_index] = bounds;
    return node_index;
}

double ContactSearchTree::distance_squared(std::size_t node, const std::array<double, 3>& point) const {
    double result = 0.0;
    for (std::size_t component = 0; component < 3; ++component) {
        const double distance =
            point[component] < _nodes[node].minimum[component]   ? _nodes[node].minimum[component] - point[component]
            : point[component] > _nodes[node].maximum[component] ? point[component] - _nodes[node].maximum[component]
                                                                 : 0.0;
        result += distance * distance;
    }
    return result;
}

void ContactSearchTree::push_node(ContactSearchQuery& query, std::size_t node) const {
    query.queue.push_back({node, distance_squared(node, query.point)});
    std::push_heap(query.queue.begin(), query.queue.end(), farther_contact_search_entry);
}

void ContactSearchTree::begin_query(const std::array<double, 3>& point, ContactSearchQuery& query) const {
    query.point = point;
    query.queue.clear();
    if (!_nodes.empty()) push_node(query, 0);
}

bool ContactSearchTree::next_candidate(ContactSearchQuery& query, double maximum_distance, std::size_t& item) const {
    const double maximum_distance_squared = maximum_distance * maximum_distance;
    while (!query.queue.empty()) {
        std::pop_heap(query.queue.begin(), query.queue.end(), farther_contact_search_entry);
        const ContactSearchQueueEntry entry = query.queue.back();
        query.queue.pop_back();
        if (entry.distance_squared > maximum_distance_squared) {
            query.queue.clear();
            return false;
        }
        const Node& node = _nodes[entry.node];
        if (node.item != std::numeric_limits<std::size_t>::max()) {
            item = node.item;
            return true;
        }
        push_node(query, node.left);
        push_node(query, node.right);
    }
    return false;
}
} // namespace fuelsim::spatial_detail

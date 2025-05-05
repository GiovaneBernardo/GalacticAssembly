#pragma once
#include "pch.h"

namespace Plaza {
	struct Vec3Comparator {
		bool operator()(const glm::vec3& a, const glm::vec3& b) const {
			if (a.x != b.x) return a.x < b.x;
			if (a.y != b.y) return a.y < b.y;
			return a.z < b.z;
		}
	};
	struct IVec3Comparator {
		std::size_t operator()(const glm::ivec3& key) const noexcept {
			size_t hx = std::hash<int>{}(key.x);
			size_t hy = std::hash<int>{}(key.y);
			size_t hz = std::hash<int>{}(key.z);
			return hx ^ (hy << 1) ^ (hz << 2);
		}
	};
}
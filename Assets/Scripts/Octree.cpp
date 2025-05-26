#include "Octree.h"

namespace Plaza {
	bool RayIntersectsAABB(const glm::vec3& rayOrigin, const glm::vec3& rayDir, const glm::vec3& aabbMin, const glm::vec3& aabbMax, float& tmin) {
		tmin = 0.0f;
		float tmax = std::numeric_limits<float>::max();

		for (int i = 0; i < 3; i++) {
			if (fabs(rayDir[i]) < 1e-6f) {
				if (rayOrigin[i] < aabbMin[i] || rayOrigin[i] > aabbMax[i])
					return false;
			} else {
				float ood = 1.0f / rayDir[i];
				float t1 = (aabbMin[i] - rayOrigin[i]) * ood;
				float t2 = (aabbMax[i] - rayOrigin[i]) * ood;
				if (t1 > t2) std::swap(t1, t2);
				tmin = glm::max(tmin, t1);
				tmax = glm::min(tmax, t2);
				if (tmin > tmax) return false;
			}
		}
		return true;
	}

	OctreeNode* OctreeNode::TraverseRay(OctreeNode* node, const glm::vec3& rayOrigin, const glm::vec3& rayDirection) {
		if (!node) return nullptr;

		glm::vec3 min = node->origin - node->size * 0.5f;
		glm::vec3 max = node->origin + node->size * 0.5f;
		float tmin;
		if (!RayIntersectsAABB(rayOrigin, rayDirection, min, max, tmin))
			return nullptr;

		if (node->isLeaf)
			return node;

		struct ChildHit {
			OctreeNode* child;
			float tmin;
		};

		std::vector<ChildHit> hits;
		for (const auto& childPtr : node->children) {
			if (!childPtr) continue;
			OctreeNode* child = childPtr.get();

			glm::vec3 childMin = child->origin - child->size * 0.5f;
			glm::vec3 childMax = child->origin + child->size * 0.5f;
			float childTmin;
			if (RayIntersectsAABB(rayOrigin, rayDirection, childMin, childMax, childTmin)) {
				hits.push_back({child, childTmin});
			}
		}

		std::sort(hits.begin(), hits.end(), [](const ChildHit& a, const ChildHit& b) {
			return a.tmin < b.tmin;
		});

		for (const auto& hit : hits) {
			OctreeNode* result = TraverseRay(hit.child, rayOrigin, rayDirection);
			if (result)
				return result;
		}

		return nullptr;
	}

}
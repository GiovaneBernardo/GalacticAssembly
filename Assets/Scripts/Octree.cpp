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

		for (const auto& child : node->children) {
			if (child) {
				OctreeNode* result = TraverseRay(child.get(), rayOrigin, rayDirection);
				if (result)
					return result;
			}
		}

		return nullptr;
	}

}
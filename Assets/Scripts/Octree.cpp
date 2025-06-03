#include "Octree.h"

#include "Constants.h"

namespace Plaza {
	bool RayIntersectsAABB(const glm::vec3& rayOrigin, const glm::vec3& rayDir, const glm::vec3& aabbMin,
						   const glm::vec3& aabbMax, float& tmin) {
		tmin = 0.0f;
		float tmax = std::numeric_limits<float>::max();

		for (int i = 0; i < 3; i++) {
			if (fabs(rayDir[i]) < 1e-6f) {
				if (rayOrigin[i] < aabbMin[i] || rayOrigin[i] > aabbMax[i])
					return false;
			}
			else {
				float ood = 1.0f / rayDir[i];
				float t1 = (aabbMin[i] - rayOrigin[i]) * ood;
				float t2 = (aabbMax[i] - rayOrigin[i]) * ood;
				if (t1 > t2)
					std::swap(t1, t2);
				tmin = glm::max(tmin, t1);
				tmax = glm::min(tmax, t2);
				if (tmin > tmax)
					return false;
			}
		}
		return true;
	}

	OctreeNode* OctreeNode::TraverseRay(OctreeNode* node, const glm::vec3& rayOrigin, const glm::vec3& rayDirection) {
		if (!node)
			return nullptr;

		glm::vec3 min = glm::vec3(node->origin) - node->size * 0.5f;
		glm::vec3 max = glm::vec3(node->origin) + node->size * 0.5f;
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
			if (!childPtr)
				continue;
			OctreeNode* child = childPtr.get();

			glm::vec3 childMin = glm::vec3(child->origin) - child->size * 0.5f;
			glm::vec3 childMax = glm::vec3(child->origin) + child->size * 0.5f;
			float childTmin;
			if (RayIntersectsAABB(rayOrigin, rayDirection, childMin, childMax, childTmin)) {
				hits.push_back({child, childTmin});
			}
		}

		std::sort(hits.begin(), hits.end(), [](const ChildHit& a, const ChildHit& b) { return a.tmin < b.tmin; });

		for (const auto& hit : hits) {
			OctreeNode* result = TraverseRay(hit.child, rayOrigin, rayDirection);
			if (result)
				return result;
		}

		return nullptr;
	}

	const glm::vec3 offsets[8] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
								  {-1, -1, 1},	{1, -1, 1},	 {-1, 1, 1},  {1, 1, 1}};

	void OctreeNode::SubDivide() {
		float childSize = size / 2.0f;

		for (int i = 0; i < 8; ++i) {
			glm::vec3 offset = offsets[i];
			glm::vec3 childOrigin = glm::vec3(origin) + offset * (childSize / 2.0f);

			// Snap to voxel grid here:
			float voxelSize = childSize / SCALAR_FIELD_SIZE;
			childOrigin = glm::round(childOrigin / voxelSize) * voxelSize;

			children[i] = std::make_unique<OctreeNode>();
			children[i]->origin = childOrigin;
			children[i]->size = childSize;
			children[i]->isLeaf = true;
		}

		isLeaf = false;
	}

	void OctreeNode::Update(const glm::vec3& playerPosition, const float planetRadius, const bool generateSmallestLOD) {
		// if (IsCompletelyOutsideSphere(glm::vec3(0.0f), planetRadius)) {
		//	// Don't subdivide or update
		//	PL_INFO("Outside sphere");
		//	return;
		// }
		if (glm::length(glm::vec3(origin)) - size * 0.5f > planetRadius)
			return;

		glm::vec3 halfSize = glm::vec3(size / 2.0f);
		glm::vec3 minBounds = glm::vec3(origin) - halfSize;
		glm::vec3 maxBounds = glm::vec3(origin) + halfSize;

		bool inside = playerPosition.x >= minBounds.x && playerPosition.x <= maxBounds.x &&
					  playerPosition.y >= minBounds.y && playerPosition.y <= maxBounds.y &&
					  playerPosition.z >= minBounds.z && playerPosition.z <= maxBounds.z;

		if (size == planetRadius)
			inside = true;

		float minSize = SCALAR_FIELD_SIZE;
		if ((inside && size > minSize) || (generateSmallestLOD && size > minSize)) {
			if (isLeaf)
				SubDivide();
			for (auto& child : children) {
				if (child)
					child->Update(playerPosition, planetRadius, generateSmallestLOD);
			}
			isLeaf = false;
		}
		else if (!inside || size <= minSize) {
			// Too small for current distance — collapse node
			for (auto& child : children) {
				child.reset();
			}
			isLeaf = true;
		}
	}

	bool OctreeNode::IsCompletelyOutsideSphere(const glm::vec3& center, float radius) {
		glm::vec3 halfSize = glm::vec3(size / 2.0f);
		glm::vec3 minBounds = glm::vec3(origin) - halfSize;
		glm::vec3 maxBounds = glm::vec3(origin) + halfSize;

		// Clamp center to AABB
		glm::vec3 closestPoint = glm::clamp(center, minBounds, maxBounds);
		float distance = glm::distance(center, closestPoint);

		return distance > radius;
	}
} // namespace Plaza
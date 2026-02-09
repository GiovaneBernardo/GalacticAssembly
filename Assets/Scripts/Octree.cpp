#include "Octree.h"
#include "pch.h"

#include "Constants.h"

namespace Plaza {
	const glm::vec3 offsets[8] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
								  {-1, -1, 1},	{1, -1, 1},	 {-1, 1, 1},  {1, 1, 1}};

	void PlanetOctreeNode::SubDivide() {
		float childSize = mSize / 2.0f;

		for (int i = 0; i < 8; ++i) {
			glm::vec3 offset = offsets[i];
			glm::vec3 childOrigin = glm::vec3(mOrigin) + offset * (childSize / 2.0f);

			// Snap to voxel grid here:
			float voxelSize = childSize / SCALAR_FIELD_SIZE;
			childOrigin = glm::round(childOrigin / voxelSize) * voxelSize;

			mChildren[i] = std::make_unique<PlanetOctreeNode>();
			mChildren[i]->mOrigin = childOrigin;
			mChildren[i]->mSize = childSize;
			mChildren[i]->mIsLeaf = true;
		}

		mIsLeaf = false;
	}

	void PlanetOctreeNode::Update(const glm::vec3& playerPosition, const float planetRadius,
								  const bool generateSmallestLOD) {
		if (glm::length(glm::vec3(mOrigin)) - mSize * 0.5f > planetRadius)
			return;

		glm::vec3 halfSize = glm::vec3(mSize / 2.0f);
		glm::vec3 minBounds = glm::vec3(mOrigin) - halfSize;
		glm::vec3 maxBounds = glm::vec3(mOrigin) + halfSize;

		bool inside = playerPosition.x >= minBounds.x && playerPosition.x <= maxBounds.x &&
					  playerPosition.y >= minBounds.y && playerPosition.y <= maxBounds.y &&
					  playerPosition.z >= minBounds.z && playerPosition.z <= maxBounds.z;

		if (mSize == planetRadius)
			inside = true;

		float minSize = SCALAR_FIELD_SIZE;
		if ((inside && mSize > minSize) || (generateSmallestLOD && mSize > minSize)) {
			if (mIsLeaf)
				SubDivide();
			for (auto& child : mChildren) {
				if (child)
					child->Update(playerPosition, planetRadius, generateSmallestLOD);
			}
			mIsLeaf = false;
		}
		else if (!inside || mSize <= minSize) {
			// Too small for current distance — collapse node
			for (auto& child : mChildren) {
				child.reset();
			}
			mIsLeaf = true;
		}
	}

	void PlanetOctreeNode::GetLeafNodes(std::vector<PlanetOctreeNode*>& nodes) {
		if (mIsLeaf)
			nodes.push_back(this);
		else {
			for (auto& child : mChildren) {
				if (child) {
					child->GetLeafNodes(nodes);
				}
			}
		}
	}
} // namespace Plaza
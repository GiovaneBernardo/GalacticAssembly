#pragma once
#include "pch.h"
#include "Utils/Comparators.h"

namespace Plaza {
	typedef std::unordered_map<glm::ivec3, float, IVec3Comparator> ScalarField;

	class PlanetOctreeNode {
	  public:
		glm::ivec3 mOrigin;
		int mSize;
		std::array<std::unique_ptr<PlanetOctreeNode>, 8> mChildren;
		bool mIsLeaf = true;
		ScalarField mScalarField;
		uint64_t mEntityUuid;

		virtual ~PlanetOctreeNode() = default;
		PlanetOctreeNode() = default;
		PlanetOctreeNode(PlanetOctreeNode&&) = default;
		PlanetOctreeNode& operator=(PlanetOctreeNode&&) = default;
		virtual void SubDivide();
		void Update(const glm::vec3& playerPosition, const float planetRadius, const bool generateSmallestLOD);
		bool IsCompletelyOutsideSphere(const glm::vec3& center, float radius);
		void GetLeafNodes(std::vector<PlanetOctreeNode*>& nodes);

		static PlanetOctreeNode* TraverseRay(PlanetOctreeNode* node, const glm::vec3& rayOrigin,
											 const glm::vec3& rayDirection);
	};
} // namespace Plaza
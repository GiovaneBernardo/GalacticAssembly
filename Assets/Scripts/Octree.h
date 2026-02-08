#pragma once
#include "pch.h"
#include "Utils/Comparators.h"

namespace Plaza {
    typedef std::unordered_map<glm::ivec3, float, IVec3Comparator> ScalarField;

    class OctreeNode {
    public:
        glm::ivec3 origin;
        int size;
        std::array<std::unique_ptr<OctreeNode>, 8> children;
        bool isLeaf = true;
        ScalarField scalarField;
    	uint64_t entityUuid;

        virtual ~OctreeNode() = default;
        OctreeNode() = default;
        OctreeNode(OctreeNode&&) = default;
        OctreeNode& operator=(OctreeNode&&) = default;
        virtual void SubDivide();
        void Update(const glm::vec3& playerPosition, const float planetRadius, const bool generateSmallestLOD);
        bool IsCompletelyOutsideSphere(const glm::vec3& center, float radius);

    	static OctreeNode* TraverseRay(OctreeNode* node, const glm::vec3& rayOrigin, const glm::vec3& rayDirection);
    };
}
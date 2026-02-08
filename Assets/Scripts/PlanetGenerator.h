#pragma once
#include "pch.h"
#include "Octree.h"
#include "Utils/Comparators.h"

namespace Plaza {

	//constexpr float PLANET_RADIUS = 1000.0f;
	//constexpr float ROOT_SIZE = PLANET_RADIUS * 2.0f; // Full diameter
	constexpr int FIELD_RES = 32;
	constexpr int MIN_NODE_SIZE = 32;
	//constexpr float ROOT_VOXEL_SIZE = ROOT_SIZE / FIELD_RES;

	// Forward declarations
	enum class TransvoxelFace : int;
	struct EnhancedOctreeNode;

	class PlanetGenerator : public CppScript {
	public:
		int mGridSize = 64;
		float mIsoLevel = 0.0f;
		float mNoiseFrequency = 0.0005f;
		float mNoiseAmplitude = 0.1f;
		float mMaxHeight = 250.0f;
		float mMinHeight = 0.0f;
		int mScale = 1.0f;
		glm::vec3 mCubicChunksSize = glm::vec3(256.0f, 256.0f, 256.0f);
		int mPlanetRadius = 1024;
		bool mGenerateSmallestLOD;
		std::vector<Material*> mMaterialsVector;

		PerlinNoise mPerlin;

		void OnStart(Scene* scene);
		void OnUpdateEditorGUI(Scene* scene);

		std::unordered_map<glm::vec3, float, Vec3Comparator> mLastGeneratedGrid;
		Entity* mLastGeneratedPlanet = nullptr;

		Entity* mPlanetEntity = nullptr;

		ScalarField GenerateSphere(const glm::vec3& origin, float size, float planetRadius, PerlinNoise& perlin);

		void SubdivideChunk(Scene* scene, const glm::ivec3& offset, int size, int depth, int maxDepth, PerlinNoise& perlin);
		void GeneratePlanetChunks(Scene* scene, const int chunkSize, const int gridSize, const float isoLevel, const int scale);

		static void UpdateNodeMesh(Scene* scene, OctreeNode* node);

		static std::unordered_map<glm::ivec3, float, IVec3Comparator> GenerateChunkGrid(PerlinNoise& perlin, const glm::ivec3& chunkOffset, const int resolution, const int planetRadius, const int chunkSize, const float isoLevel, const float noiseFrequency, const float noiseAmplitude);
	private:
		Entity* GeneratePlanet(Scene* scene);
		Entity* GeneratePlanet(Scene* scene, const std::unordered_map<glm::ivec3, float, IVec3Comparator>& grid, const int gridSize, const float isoLevel);
		static Mesh* GenerateMesh(const std::unordered_map<glm::ivec3, float, IVec3Comparator>& grid, const glm::vec3& offset, const int gridSize, const float isoLevel, const float resolution, const glm::ivec3 chunkOffset);
		std::unordered_map<glm::ivec3, float, IVec3Comparator> GenerateGrid(const int gridSize, const float isoLevel, const float noiseFrequency, const float noiseAmplitude);

		void BuildNeighborRelationships(EnhancedOctreeNode& root);
		void CollectAllNodes(EnhancedOctreeNode& node, std::vector<EnhancedOctreeNode*>& outNodes);
		EnhancedOctreeNode* FindNeighborAtFace(EnhancedOctreeNode* node, TransvoxelFace face,
															const std::vector<EnhancedOctreeNode*>& allNodes);
		glm::vec3 GetFaceOffset(TransvoxelFace face);
		void GenerateInteriorMesh(const ScalarField& grid, EnhancedOctreeNode* node, float isoLevel,
											  float resolution, std::vector<glm::vec3>& vertices,
											  std::vector<unsigned int>& indices, std::vector<glm::vec3>& normals,
											  std::vector<glm::vec2>& uvs, std::vector<unsigned int>& materials);
		void GenerateTransitionCells(const ScalarField& highResField, const ScalarField& lowResField,
												 TransvoxelFace face, EnhancedOctreeNode* highResNode,
												 EnhancedOctreeNode* lowResNode, float isoLevel, float resolution,
												 std::vector<glm::vec3>& vertices, std::vector<unsigned int>& indices,
												 std::vector<glm::vec3>& normals, std::vector<glm::vec2>& uvs,
												 std::vector<unsigned int>& materials);
		void GetFaceCoordinates(TransvoxelFace face, int u, int v, glm::ivec3 coords[4]);
		void GenerateTransitionQuad(float samples[4], glm::ivec3 coords[4], TransvoxelFace face,
												EnhancedOctreeNode* node, float isoLevel, float resolution,
												std::vector<glm::vec3>& vertices, std::vector<unsigned int>& indices,
												std::vector<glm::vec3>& normals, std::vector<glm::vec2>& uvs,
												std::vector<unsigned int>& materials);
		Mesh* GenerateTransvoxelMesh(EnhancedOctreeNode* node, float isoLevel);
	};
	PL_REGISTER_SCRIPT(PlanetGenerator);
}
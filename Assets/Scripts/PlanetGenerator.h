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

	class PlanetGenerator : public CppScript {
	public:
		int mGridSize = 64;
		float mIsoLevel = 0.8f;
		float mNoiseFrequency = 2.0f;
		float mNoiseAmplitude = 2.7f;
		float mMaxHeight = 15000.0f;
		float mMinHeight = 0.0f;
		int mScale = 1.0f;
		glm::vec3 mCubicChunksSize = glm::vec3(256.0f, 256.0f, 256.0f);
		float mPlanetRadius = 1000.0f;

		PerlinNoise mPerlin;

		void OnStart(Scene* scene);
		void OnUpdateEditorGUI(Scene* scene);

		std::unordered_map<glm::vec3, float, Vec3Comparator> mLastGeneratedGrid;
		Entity* mLastGeneratedPlanet = nullptr;

		Entity* mPlanetEntity = nullptr;

		ScalarField GenerateSphere(const glm::vec3& origin, float size, float planetRadius, PerlinNoise& perlin);

		void SubdivideChunk(Scene* scene, const glm::ivec3& offset, int size, int depth, int maxDepth, PerlinNoise& perlin);
		void GeneratePlanetChunks(Scene* scene, const int chunkSize, const int gridSize, const float isoLevel, const int scale);

		static std::unordered_map<glm::ivec3, float, IVec3Comparator> GenerateChunkGrid(PerlinNoise& perlin, const glm::ivec3& chunkOffset, const int resolution, const int planetRadius, const int chunkSize, const float isoLevel, const float noiseFrequency, const float noiseAmplitude);
	private:
		Entity* GeneratePlanet(Scene* scene, const std::unordered_map<glm::ivec3, float, IVec3Comparator>& grid, const int gridSize, const float isoLevel);
		Mesh* GenerateMesh(const std::unordered_map<glm::ivec3, float, IVec3Comparator>& grid, const glm::vec3& offset, const int gridSize, const float isoLevel, const float resolution, const glm::ivec3 chunkOffset);
		std::unordered_map<glm::ivec3, float, IVec3Comparator> GenerateGrid(const int gridSize, const float isoLevel, const float noiseFrequency, const float noiseAmplitude);
	};
	PL_REGISTER_SCRIPT(PlanetGenerator);
}
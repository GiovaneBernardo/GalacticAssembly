#pragma once
#include "pch.h"
#include "Octree.h"
#include "Utils/Comparators.h"

namespace Plaza {

	constexpr int FIELD_RES = 32;
	constexpr int MIN_NODE_SIZE = 32;

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

	private:
		Entity* GeneratePlanet(Scene* scene);
	};
	PL_REGISTER_SCRIPT(PlanetGenerator);
}

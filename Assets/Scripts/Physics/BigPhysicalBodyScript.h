#pragma once
#include "pch.h"
#include "Utils/Comparators.h"

namespace Plaza {
	class BigPhysicalBodyScript : public CppScript {
	public:
		void OnStart(Scene* scene);
		void OnUpdate(Scene* scene);

		glm::vec3 mLinearVelocity = glm::vec3(0.0f);
		glm::vec3 mAngularVelocity = glm::vec3(0.0f);

		struct Chunk {
			std::unordered_map<glm::ivec3, float, IVec3Comparator> mScalarField;
		};

		std::unordered_map<glm::ivec3, Chunk, IVec3Comparator> mChunks;

		float mRadius;
		float mWeight;
		glm::vec3 mPosition;
		glm::quat mRotation;
		float maxResolution = 1.0f;
		float minResolution = 8.0f;

		float mGravitationalAcceleration = 9.81f;
	private:
	};
	PL_REGISTER_SCRIPT(BigPhysicalBodyScript);
}
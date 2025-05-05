#pragma once
#include "pch.h"
#include "BigPhysicalBodyScript.h"

namespace Plaza {
	class PhysicalBodyScript : public CppScript {
	public:
		void OnStart(Scene* scene);
		void OnUpdate(Scene* scene);

		glm::vec3 mLinearVelocity = glm::vec3(0.0f);
		glm::vec3 mAngularVelocity = glm::vec3(0.0f);

		static inline std::map<uint64_t, BigPhysicalBodyScript*> mBigPhysicalBodies = std::map<uint64_t, BigPhysicalBodyScript*>();
	private:
	};
	PL_REGISTER_SCRIPT(PhysicalBodyScript);
}
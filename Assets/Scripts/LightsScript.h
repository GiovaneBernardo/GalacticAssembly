#pragma once
#include "pch.h"

namespace Plaza {
	class LightsScript : public CppScript {
	public:
		struct LightTransform {
			TransformComponent* mTransform = nullptr;
			glm::vec3 mMovingDirection = glm::vec3(1.0f, 0.0f, 0.0f);
			float mSpeed = 1.0f;
		};
		std::vector<LightTransform> mLights = std::vector<LightTransform>();

		glm::vec3 mOrigin = glm::vec3(0.0f);
		glm::vec3 mPoint0 = glm::vec3(0.0f);
		glm::vec3 mPoint1 = glm::vec3(0.0f);
		void OnStart(Scene* scene);
		void OnUpdate(Scene* scene);

	private:
	};
	PL_REGISTER_SCRIPT(LightsScript);
}
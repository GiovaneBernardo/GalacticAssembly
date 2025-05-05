#pragma once
#include "pch.h"

namespace Plaza {
	class CameraScript : public CppScript {
	public:
		glm::vec3 mVelocity = glm::vec3(0.0f);
		std::vector<glm::mat4> mMatrices = std::vector<glm::mat4>();
		int mPosCount = 9;
		float mSpeed = 10.0f;
		float mCurrentTime = 0.0f;
		int mCurrentPosIndex = 0;
		void OnStart(Scene* scene);
		void OnUpdate(Scene* scene);

	private:
	};
	PL_REGISTER_SCRIPT(CameraScript);
}
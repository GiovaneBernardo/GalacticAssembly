#pragma once
#include "pch.h"

namespace Plaza {
	class CameraScript : public CppScript {
	public:
		glm::vec3 mVelocity = glm::vec3(0.0f);
		float mMaxDeceleration = 8.0f;
		std::vector<glm::mat4> mMatrices = std::vector<glm::mat4>();
		int mPosCount = 9;
		float mSpeed = 10.0f;
		float mCurrentTime = 0.0f;
		int mCurrentPosIndex = 0;
		bool mIsFlying = false;
		bool mDampers = true;
		void OnStart(Scene* scene);
		void OnUpdate(Scene* scene);

	private:
		void ModifyTerrain(Scene* scene);
	};
	PL_REGISTER_SCRIPT(CameraScript);
}
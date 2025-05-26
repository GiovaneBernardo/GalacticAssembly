#pragma once
#include "pch.h"

namespace Plaza {
	class CameraScript : public CppScript {
	public:
		glm::vec3 mVelocity = glm::vec3(0.0f);
		float mMaxDeceleration = 512.0f;
		float mAcceleration = 32.0f;
		bool mShowAffectedNode = false;
		bool mShowDebugEntities = false;
		std::vector<glm::mat4> mMatrices = std::vector<glm::mat4>();
		Material* mDebugMaterial;
		int mPosCount = 9;
		float mSpeed = 10.0f;
		float mCurrentTime = 0.0f;
		int mCurrentPosIndex = 0;
		bool mIsFlying = false;
		bool mDampers = true;
		float mDigRadius = 7.0f;
		float mDigStrength = 1.0f;
		void OnStart(Scene* scene);
		void OnUpdate(Scene* scene);
		void OnUpdateEditorGUI(Scene* scene);

	private:
		void ModifyTerrain(Scene* scene);
	};
	PL_REGISTER_SCRIPT(CameraScript);
}
#include "LightsScript.h"

namespace Plaza {
	void LightsScript::OnStart(Scene* scene) {
		mOrigin = glm::vec3(0.0f, 0.0f, -88.0f);
		mPoint0 = mOrigin - glm::vec3(13, 0.0f, 5.7f);
		mPoint1 = mOrigin + glm::vec3(13, 10.0f, 5.7f);

		static  std::random_device rd2;
		static  std::mt19937 gen2(rd2());
		static  std::uniform_real_distribution<> dis2(-1.0, 1.0);
		static  std::uniform_real_distribution<> dis3(0.0, 1.0);
		static  std::uniform_real_distribution<> dis4(3.0, 10.0);
		static  std::uniform_real_distribution<> dis5(0.0, 5.0);

		for (int i = 0; i < 128; ++i) {
			Entity* newLight = Scene::GetActiveScene()->GetEntity(ECS::EntitySystem::Instantiate(scene, CppHelper::FindEntity("LightInstance")->uuid));
			mLights.push_back(LightTransform(scene->GetComponent<TransformComponent>(newLight->uuid), glm::vec3(dis2(gen2), dis2(gen2), dis2(gen2)), dis5(gen2) + 1));
			scene->GetComponent<Light>(newLight->uuid)->color = glm::vec3(dis3(gen2), dis3(gen2), dis3(gen2));
			scene->GetComponent<Light>(newLight->uuid)->radius = dis4(gen2);
			scene->GetComponent<Light>(newLight->uuid)->intensity = dis5(gen2);
		}
	}

	bool IsInside(const glm::vec3& point0, const glm::vec3& point1, const glm::vec3& vec) {
		return vec.x > point0.x && vec.x < point1.x && vec.y > point0.y && vec.y < point1.y && vec.z > point0.z && vec.z < point1.z;
	}

	void LightsScript::OnUpdate(Scene* scene) {
		static  std::random_device rd;
		static  std::mt19937 gen(rd());
		static  std::uniform_real_distribution<> dis(-0.3, 0.3);
		for (LightTransform& lightTransform : mLights) {
			glm::vec3 position = lightTransform.mTransform->GetWorldPosition();
			if (position.x < mPoint0.x || position.x > mPoint1.x)
				lightTransform.mMovingDirection.x *= -1.0f;
			if (position.y < mPoint0.y || position.y > mPoint1.y)
				lightTransform.mMovingDirection.y *= -1.0f;
			if (position.z < mPoint0.z || position.z > mPoint1.z)
				lightTransform.mMovingDirection.z *= -1.0f;
			//if (!IsInside(mPoint0, mPoint1, lightTransform.mTransform->GetWorldPosition())) {
			//	lightTransform.mMovingDirection = glm::normalize(-lightTransform.mMovingDirection + glm::vec3(dis(gen), dis(gen), dis(gen)));
			//}

			ECS::TransformSystem::SetLocalPosition(*lightTransform.mTransform, scene, lightTransform.mTransform->GetWorldPosition() + lightTransform.mMovingDirection * Time::GetDeltaTime() * lightTransform.mSpeed);
		}
	}
}
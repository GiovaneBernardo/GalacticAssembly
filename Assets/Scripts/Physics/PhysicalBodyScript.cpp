#include "PhysicalBodyScript.h"

namespace Plaza {
//	std::mutex PhysicalBodyScript::sMapMutex = std::mutex();

	void PhysicalBodyScript::OnStart(Scene* scene) {

	}

	void PhysicalBodyScript::OnUpdate(Scene* scene) {
		// Iterate over every big physical body and add a velocity pushing this physical body to the big physical body
		for (const auto& [key, body] : sBigPhysicalBodies) {
			glm::vec3 direction = scene->GetComponent<TransformComponent>(key)->GetWorldPosition() - scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldPosition();
			float distance = glm::length(direction);
			if (distance < 0.0001f)
				continue;
			distance = glm::max(distance, 5.0f);
			float acceleration = (6.67430e-11 * body->mWeight) / (distance * distance);
			PL_CORE_INFO("Acceleration: {}", acceleration);
			// Replace 1.0f / 60.0f by delta time -------------------------------------------------------------------------------------------------------------
			this->mLinearVelocity += glm::normalize(direction) * acceleration * glm::vec3(1.0f / 60.0f);;
		}
		// Update the position with the velocity
		const glm::vec3 correctedVelocity = mLinearVelocity * glm::vec3(1.0f / 60.0f);//glm::vec3(mLinearVelocity.x * Time::deltaTime, mLinearVelocity.y * Time::deltaTime, mLinearVelocity.z * Time::deltaTime);
		ECS::TransformSystem::SetWorldPosition(*scene->GetComponent<TransformComponent>(this->mEntityUuid), scene, scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldPosition() + correctedVelocity);

		PL_CORE_INFO("X: {} \n Y: {} \n Z: {}", mLinearVelocity.x, mLinearVelocity.y, mLinearVelocity.z);
		//this->mLinearVelocity = glm::vec3(0.0f);
		// Reset the velocity if it is colliding with something
		//RaycastHit hit;
		//Physics::Raycast(scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldPosition(), glm::normalize(scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldRotation() * glm::vec3(0.0f, -1.0f, 0.0f)), hit);
		//if (!hit.missed) {
		//	PL_CORE_INFO("Hitting");
		//	mLinearVelocity = glm::vec3(0.0f);
		//}
	}


}
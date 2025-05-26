#include "PhysicalBodyScript.h"

namespace Plaza {
	//	std::mutex PhysicalBodyScript::sMapMutex = std::mutex();

	void PhysicalBodyScript::OnStart(Scene* scene) {}

	void PhysicalBodyScript::OnUpdate(Scene* scene) {
		// Iterate over every big physical body and add a velocity pushing this physical body to the big physical body
		for (const auto& [key, bigBody] : sBigPhysicalBodies) {
			glm::vec3 direction = scene->GetComponent<TransformComponent>(key)->GetWorldPosition() -
								  scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldPosition();
			double distance = glm::length(direction);
			if (distance < 0.0001f)
				continue;
			distance = glm::max(distance, 5.0);
			double softening = 1.0f;
			float acceleration =
				(6.67430e-11 * bigBody->mWeight * 50.00) / (distance * distance + softening * softening);
			// Replace 1.0f / 60.0f by delta time
			// -------------------------------------------------------------------------------------------------------------
			this->mLinearVelocity += glm::normalize(direction) * (acceleration / 50.0f) * glm::vec3(1.0f / 60.0f);
			if (scene->HasComponent<RigidBody>(this->mEntityUuid))
				scene->GetComponent<RigidBody>(this->mEntityUuid)
					->SetVelocity(scene->GetComponent<RigidBody>(this->mEntityUuid)->GetVelocity() +
								  glm::normalize(direction) * (acceleration / 50.0f) * glm::vec3(1.0f / 60.0f));
		}
		// Update the position with the velocity
		const glm::vec3 correctedVelocity =
			mLinearVelocity *
			glm::vec3(1.0f / 60.0f); // glm::vec3(mLinearVelocity.x * Time::deltaTime, mLinearVelocity.y *
									 // Time::deltaTime, mLinearVelocity.z * Time::deltaTime);
		//ECS::TransformSystem::SetWorldPosition(
		//	*scene->GetComponent<TransformComponent>(this->mEntityUuid), scene,
		//	scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldPosition() + correctedVelocity);

		// this->mLinearVelocity = glm::vec3(0.0f);
		//  Reset the velocity if it is colliding with something
		// RaycastHit hit;
		// Physics::Raycast(scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldPosition(),
		// glm::normalize(scene->GetComponent<TransformComponent>(this->mEntityUuid)->GetWorldRotation() *
		// glm::vec3(0.0f, -1.0f, 0.0f)), hit); if (!hit.missed) { 	PL_CORE_INFO("Hitting"); 	mLinearVelocity =
		// glm::vec3(0.0f);
		// }
	}

	void PhysicalBodyScript::OnCollide(Scene* scene, const RaycastHit hit) {
		mLinearVelocity = glm::vec3(0.0f);
		scene->GetComponent<RigidBody>(this->mEntityUuid)->SetVelocity(glm::vec3(0.0f, 0.0f, 0.0f));
	}

} // namespace Plaza
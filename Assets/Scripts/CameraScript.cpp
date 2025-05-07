#include "CameraScript.h"
#include "Engine/Core/Time.h"

namespace Plaza {
	void CameraScript::OnStart(Scene* scene) {

	}

    float ApplyDragToFloat(float value, const float drag) {
		const float realDrag = glm::max(drag, value * drag);
	    if (value >= 0.0f)
	        return glm::max(0.0f, value - realDrag);

	    return glm::min(0.0f, value + realDrag);
	}

    glm::vec3 ApplyDragToVec3(glm::vec3 value, const float drag) {
	    return glm::vec3(
            ApplyDragToFloat(value.x, drag),
            ApplyDragToFloat(value.y, drag),
            ApplyDragToFloat(value.z, drag));
	}

    void CameraScript::OnUpdate(Scene* scene) {
        scene = Scene::GetActiveScene();

        uint64_t playerUuid = CppHelper::FindEntity("PlayerEntity")->uuid;
        TransformComponent* playerTransform = scene->GetComponent<TransformComponent>(playerUuid);

        // Rotate Player
        glm::vec2 mouseDelta = Input::Cursor::GetDeltaMousePosition();
        float sensitivity = 0.2f * Time::GetDeltaTime();
        glm::vec3 eulerDelta(
            mouseDelta.y * sensitivity,  // Pitch
            -mouseDelta.x * sensitivity, // Yaw
            0.0f                         // Roll
        );
        glm::quat deltaRotation = glm::quat(eulerDelta);
        glm::quat currentRotation = playerTransform->GetLocalQuaternion();
        glm::quat newRotation = currentRotation * deltaRotation;

        // Move Player
		auto& transform = *scene->GetComponent<TransformComponent>(playerUuid);
		glm::quat rotation = transform.GetWorldQuaternion();
		float dt = Time::GetDeltaTime();

		if (Input::GetKeyDown(GLFW_KEY_W))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, 0.0f, 1.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_S))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, 0.0f, -1.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_D))
			mVelocity += glm::rotate(rotation, glm::vec3(-1.0f, 0.0f, 0.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_A))
			mVelocity += glm::rotate(rotation, glm::vec3(1.0f, 0.0f, 0.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_SPACE))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, 1.0f, 0.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_C))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, -1.0f, 0.0f)) * dt;

	    float drag = 0.8f * Time::GetDeltaTime();
	    //mVelocity = ApplyDragToVec3(mVelocity, drag);

        ECS::TransformSystem::SetLocalRotation(*playerTransform, scene, newRotation);
        ECS::TransformSystem::SetLocalPosition(*playerTransform, scene, playerTransform->GetWorldPosition() + mVelocity);
        ECS::TransformSystem::UpdateSelfAndChildrenTransform(*scene->GetComponent<TransformComponent>(playerUuid), nullptr, scene, true);
		mVelocity = glm::vec3(0.0f);
    }
}
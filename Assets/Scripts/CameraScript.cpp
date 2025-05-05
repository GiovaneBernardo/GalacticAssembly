#include "CameraScript.h"
#include "Engine/Core/Time.h"

namespace Plaza {
	void CameraScript::OnStart(Scene* scene) {

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
        if (Input::GetKeyDown(GLFW_KEY_W))
            mVelocity += glm::vec3(0.0f, 0.0f, 1.0f * Time::GetDeltaTime());

        ECS::TransformSystem::SetLocalRotation(*playerTransform, scene, newRotation);
        ECS::TransformSystem::SetLocalPosition(*playerTransform, scene, playerTransform->GetWorldPosition() + mVelocity);
        ECS::TransformSystem::UpdateSelfAndChildrenTransform(*scene->GetComponent<TransformComponent>(playerUuid), nullptr, scene, true);
    }
}
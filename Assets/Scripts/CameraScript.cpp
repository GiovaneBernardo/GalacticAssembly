#include "CameraScript.h"

#include "PlanetGenerator.h"
#include "Engine/Application/Callbacks/CallbacksHeader.h"
#include "Engine/Core/Time.h"
#include "Physics/BigPhysicalBodyScript.h"
#include "Constants.h"

namespace Plaza {
	void CameraScript::OnStart(Scene* scene) {
		mDebugMaterial = new Material();
		mDebugMaterial->diffuse->rgba = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
		AssetsManager::AddMaterial(mDebugMaterial);
	}

	void CameraScript::OnUpdateEditorGUI(Scene* scene) {
		ImGui::DragFloat("Dig Radius", &mDigRadius);
		ImGui::DragFloat("Dig Strength", &mDigStrength);
		ImGui::DragFloat("Acceleration", &mAcceleration);
		ImGui::DragFloat("Deceleration", &mMaxDeceleration);
		ImGui::Checkbox("Show Affected Node", &mShowAffectedNode);
		ImGui::Checkbox("Show Debug Entities", &mShowDebugEntities);
	}

	float ApplyDragToFloat(float value, const float drag) {
		const float realDrag = glm::max(drag, value * drag);
		if (value >= 0.0f)
			return glm::max(0.0f, value - realDrag);

		return glm::min(0.0f, value + realDrag);
	}

	glm::vec3 ApplyDragToVec3(glm::vec3 value, const float drag) {
		return glm::vec3(ApplyDragToFloat(value.x, drag), ApplyDragToFloat(value.y, drag),
						 ApplyDragToFloat(value.z, drag));
	}

	void CameraScript::OnUpdate(Scene* scene) {
		scene = Scene::GetActiveScene();

		uint64_t cameraUuid = CppHelper::FindEntity("CameraEntity")->uuid;
		uint64_t playerUuid = CppHelper::FindEntity("PlayerEntity")->uuid;
		TransformComponent* playerTransform = scene->GetComponent<TransformComponent>(playerUuid);

		// Move Player
		TransformComponent* cameraTransform = scene->GetComponent<TransformComponent>(cameraUuid);
		glm::quat rotation;
		rotation = playerTransform->GetWorldQuaternion();
		float dt = Time::GetDeltaTime();

		if (Input::GetKeyDown(GLFW_KEY_W))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, 0.0f, mAcceleration)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_S))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, 0.0f, -mAcceleration)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_D))
			mVelocity += glm::rotate(rotation, glm::vec3(-mAcceleration, 0.0f, 0.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_A))
			mVelocity += glm::rotate(rotation, glm::vec3(mAcceleration, 0.0f, 0.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_SPACE))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, mAcceleration, 0.0f)) * dt;
		if (Input::GetKeyDown(GLFW_KEY_C))
			mVelocity += glm::rotate(rotation, glm::vec3(0.0f, -mAcceleration, 0.0f)) * dt;

		if (mDampers && glm::distance(mVelocity, glm::vec3(0.0f)) == 0.0f) {
			const glm::vec3& velocity = scene->GetComponent<RigidBody>(this->mEntityUuid)->GetVelocity();
			mVelocity -= glm::normalize(scene->GetComponent<RigidBody>(this->mEntityUuid)->GetVelocity()) *
						 glm::smoothstep(mMaxDeceleration, 0.0f, glm::distance(velocity, glm::vec3(0.0f))) * dt;
		}

		if (scene->HasComponent<RigidBody>(this->mEntityUuid)) {
			scene->GetComponent<RigidBody>(this->mEntityUuid)
				->SetVelocity(scene->GetComponent<RigidBody>(this->mEntityUuid)->GetVelocity() + mVelocity);
		}
		else {
			ECS::TransformSystem::SetLocalPosition(*playerTransform, scene,
												   playerTransform->GetWorldPosition() + mVelocity);
		}

		// Rotate Player
		glm::vec2 mouseDelta = Input::Cursor::GetDeltaMousePosition();
		float sensitivity = 0.2f * Time::GetDeltaTime();
		glm::quat currentPlayerRotation = playerTransform->GetLocalQuaternion();

		if (mIsFlying) {
			glm::vec3 eulerDelta(mouseDelta.y * sensitivity,  // Pitch
								 -mouseDelta.x * sensitivity, // Yaw
								 0.0f						  // Roll
			);

			// Apply both X and Y to character body, camera is always looking forward here, similar to space engineers
			glm::quat deltaRotation = glm::quat(eulerDelta);
			glm::quat newPlayerRotation = currentPlayerRotation * deltaRotation;
			ECS::TransformSystem::SetLocalRotation(*playerTransform, scene, newPlayerRotation);
		}
		else {
			// Apply X to character body and Y to camera, similar to fps games
			glm::quat yawRotation = glm::angleAxis(-mouseDelta.x * sensitivity, glm::vec3(0.0f, 1.0f, 0.0f));
			glm::quat pitchRotation = glm::angleAxis(mouseDelta.y * sensitivity, glm::vec3(1.0f, 0.0f, 0.0f));

			ECS::TransformSystem::SetLocalRotation(*playerTransform, scene, yawRotation * currentPlayerRotation);

			glm::quat currentCameraRotation = cameraTransform->GetLocalQuaternion();
			ECS::TransformSystem::SetLocalRotation(*cameraTransform, scene, pitchRotation * currentCameraRotation);
		}

		float drag = 0.8f * Time::GetDeltaTime();
		// mVelocity = ApplyDragToVec3(mVelocity, drag);

		ECS::TransformSystem::UpdateSelfAndChildrenTransform(*scene->GetComponent<TransformComponent>(playerUuid),
															 nullptr, scene, true);
		mVelocity = glm::vec3(0.0f);

		if (Input::GetKeyDownOnce(GLFW_KEY_X))
			mIsFlying = !mIsFlying;
		if (Input::GetKeyDownOnce(GLFW_KEY_Z))
			mDampers = !mDampers;

		ModifyTerrain(scene);
	}

	bool ModifyScalarFieldInRadius(ScalarField& field, glm::ivec3 pos, float radius, float strength, OctreeNode* node) {
		bool modified = false;
		for (int x = -radius; x <= radius; ++x) {
			for (int y = -radius; y <= radius; ++y) {
				for (int z = -radius; z <= radius; ++z) {
					glm::ivec3 offset(x, y, z);
					glm::ivec3 targetPos = pos + offset;

					if (glm::length(glm::vec3(offset)) <= radius) {
						auto it = node->scalarField.find(targetPos);
						if (it != node->scalarField.end()) {
							float dist = glm::length(glm::vec3(offset));
							float falloff = 1.0f - (dist / radius);
							it->second += strength * falloff;
							modified = true;
						}
					}
				}
			}
		}
		return modified;
	}

	void CollectIntersectingLeafNodes(
	OctreeNode* node,
	const glm::vec3& sphereCenter,
	float sphereRadius,
	std::vector<OctreeNode*>& outNodes)
	{
		if (!node) return;

		glm::vec3 min = node->origin - glm::vec3(node->size * 0.5f);
		glm::vec3 max = node->origin + glm::vec3(node->size * 0.5f);

		// AABB-sphere intersection
		float distSq = 0.0f;
		for (int i = 0; i < 3; ++i) {
			float v = sphereCenter[i];
			if (v < min[i]) distSq += (min[i] - v) * (min[i] - v);
			if (v > max[i]) distSq += (v - max[i]) * (v - max[i]);
		}
		if (distSq > sphereRadius * sphereRadius)
			return;

		if (node->isLeaf) {
			outNodes.push_back(node);
		} else {
			for (const auto& child : node->children) {
				CollectIntersectingLeafNodes(child.get(), sphereCenter, sphereRadius, outNodes);
			}
		}
	}

	void CameraScript::ModifyTerrain(Scene* scene) {
		Entity* planet = scene->GetEntityByName("Planet");
		uint64_t cameraUuid = CppHelper::FindEntity("CameraEntity")->uuid;
		uint64_t playerUuid = CppHelper::FindEntity("PlayerEntity")->uuid;
		TransformComponent* playerTransform = scene->GetComponent<TransformComponent>(playerUuid);
		TransformComponent* cameraTransform = scene->GetComponent<TransformComponent>(cameraUuid);

		static MeshRenderer* lastMeshRenderer = nullptr;

		bool mouse0 = Input::GetMouseDown(0);
		bool mouse1 = Input::GetMouseDown(1);
		if (Input::GetMouseDown(0) && planet) {
			// Get the node the player is looking
			OctreeNode* hitNode = OctreeNode::TraverseRay(
				&static_cast<BigPhysicalBodyScript*>(scene->GetComponent<CppScriptComponent>(planet->uuid)->mScripts[0])
					 ->mRootOctreeNode,
				playerTransform->GetWorldPosition(), scene->GetComponent<Camera>(cameraTransform->mUuid)->Front);

			// Modify scalar field if node is found and it has an entity
			if (hitNode && hitNode->entityUuid) {
				// Show the node player is looking to
				if (mShowAffectedNode) {
					if (lastMeshRenderer) {
						lastMeshRenderer->ChangeMaterial(AssetsManager::GetDefaultMaterial());
					}
					lastMeshRenderer = scene->GetComponent<MeshRenderer>(hitNode->entityUuid);
					lastMeshRenderer->ChangeMaterial(mDebugMaterial);
				} else if (lastMeshRenderer) {
					lastMeshRenderer->ChangeMaterial(AssetsManager::GetDefaultMaterial());
					lastMeshRenderer = nullptr;
				}

				// Get middle of the screen
				glm::vec2 screenPos = glm::vec2(Application::Get()->appSizes->sceneSize.x * 0.5f,
												Application::Get()->appSizes->sceneSize.y * 0.5f);

				// Get depth
				float depth = VulkanRenderer::GetRenderer()
								  ->mRenderGraph->GetTexture<VulkanTexture>("SceneDepth")
								  ->ReadTexture(screenPos, sizeof(float) + sizeof(uint8_t), 1,
												VK_IMAGE_ASPECT_DEPTH_BIT | VK_IMAGE_ASPECT_STENCIL_BIT, true)
								  .x;

				// Reconstruct world position from depth
				glm::vec3 worldPosition = Renderer::ReconstructWorldPositionFromDepth(
					screenPos, Application::Get()->appSizes->sceneSize, depth, Application::Get()->activeCamera);

				if (mShowDebugEntities)
				ECS::TransformSystem::SetLocalPosition(
					*scene->GetComponent<TransformComponent>(scene->GetEntityByName("DebugNoCollider")->uuid), scene,
					worldPosition);



				//glm::vec3 local = worldPosition - hitNode->origin;
				//glm::ivec3 scalarFieldPos = glm::ivec3(local * glm::vec3(32.0f / hitNode->size));
				//glm::ivec3 scalarFieldPos = glm::ivec3(worldPosition - hitNode->origin + glm::vec3(hitNode->size / 2.0f));
				constexpr int scalarResolution = SCALAR_FIELD_SIZE;
				glm::vec3 corner = hitNode->origin - glm::vec3(hitNode->size * 0.5f);
				float voxelSize = hitNode->size / float(scalarResolution);

				if (mShowDebugEntities) {
					// Convert world position to scalar field index
					glm::vec3 local = (worldPosition - corner) / voxelSize;
					glm::ivec3 scalarFieldPos = glm::clamp(glm::ivec3(local), glm::ivec3(0), glm::ivec3(scalarResolution - 1));

					glm::vec3 debugWorldPos = corner + (glm::vec3(scalarFieldPos) + 0.5f) * voxelSize;
					ECS::TransformSystem::SetLocalPosition(
						*scene->GetComponent<TransformComponent>(scene->GetEntityByName("Sphere")->uuid),
						scene,
						debugWorldPos
					);
				}

				float strength = Input::GetMouseDown(1) ? -mDigStrength : mDigStrength;

				// Modify the scalar field points inside radius
				std::vector<OctreeNode*> affectedNodes;
				CollectIntersectingLeafNodes(
					&static_cast<BigPhysicalBodyScript*>(scene->GetComponent<CppScriptComponent>(planet->uuid)->mScripts[0])->mRootOctreeNode,
					worldPosition,
					mDigRadius,
					affectedNodes);

				for (OctreeNode* node : affectedNodes) {
					// Transform worldPos into this node scalar field space
					glm::vec3 corner = node->origin - glm::vec3(node->size * 0.5f);
					float voxelSize = node->size / float(SCALAR_FIELD_SIZE);
					glm::vec3 local = (worldPosition - corner) / voxelSize;
					glm::ivec3 scalarFieldPos = glm::clamp(glm::ivec3(local), glm::ivec3(0), glm::ivec3(SCALAR_FIELD_SIZE - 1));

					bool modified = ModifyScalarFieldInRadius(node->scalarField, scalarFieldPos, mDigRadius, strength, node);
					if (modified) {
						PlanetGenerator::UpdateNodeMesh(scene, node);
					}
				}
			}
		}
	}
} // namespace Plaza
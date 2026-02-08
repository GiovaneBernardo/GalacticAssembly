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
		glm::quat rotation = playerTransform->GetWorldQuaternion();
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
			mVelocity -= glm::normalize(velocity) *
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
		float sensitivity = 0.001f; // 0.2f * Time::GetDeltaTime();
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

	OctreeNode* FindLeafNodeAtPosition(OctreeNode* node, const glm::vec3& targetPos, float expectedSize, float epsilon = 0.001f) {
		if (!node) return nullptr;

		glm::vec3 halfSizeVec = glm::vec3(node->size * 0.5f);
		glm::vec3 min = glm::vec3(node->origin) - halfSizeVec;
		glm::vec3 max = glm::vec3(node->origin) + halfSizeVec;

		bool inside =
			targetPos.x >= min.x - epsilon && targetPos.x <= max.x + epsilon &&
			targetPos.y >= min.y - epsilon && targetPos.y <= max.y + epsilon &&
			targetPos.z >= min.z - epsilon && targetPos.z <= max.z + epsilon;

		if (!inside) return nullptr;

		if (node->isLeaf) {
			// Make sure it’s the same resolution
			if (std::abs(node->size - expectedSize) < epsilon)
				return node;
			else
				return nullptr;
		}

		for (auto& child : node->children) {
			if (child) {
				OctreeNode* found = FindLeafNodeAtPosition(child.get(), targetPos, expectedSize, epsilon);
				if (found) return found;
			}
		}

		return nullptr;
	}

	std::vector<OctreeNode*> GetNeighborNodes(OctreeNode* current, OctreeNode* root) {
		std::vector<OctreeNode*> neighbors;

		if (!current || !root) return neighbors;

		const float epsilon = 0.001f;
		const glm::vec3 origin = current->origin;
		const float size = current->size;

		// Directions to check neighbors in 6 faces
		const glm::vec3 directions[] = {
			{ size, 0.0f, 0.0f },  // +X
			{-size, 0.0f, 0.0f },  // -X
			{ 0.0f, size, 0.0f },  // +Y
			{ 0.0f,-size, 0.0f },  // -Y
			{ 0.0f, 0.0f, size },  // +Z
			{ 0.0f, 0.0f,-size }   // -Z
		};

		for (const glm::vec3& dir : directions) {
			glm::vec3 neighborCenter = origin + dir;

			// Recursively find the leaf node at this position
			OctreeNode* found = FindLeafNodeAtPosition(root, neighborCenter, size, epsilon);
			if (found && found != current) {
				neighbors.push_back(found);
			}
		}

		return neighbors;
	}

	OctreeNode* GetNodeAtOffset(OctreeNode* node, const glm::ivec3& offsetDir, OctreeNode* root) {
		if (!node || !root || glm::length(glm::vec3(offsetDir)) != 1.0f)
			return nullptr; // Only support direct axis neighbors

		constexpr float EPSILON = 0.01f;

		// Compute neighbor origin
		glm::vec3 neighborOrigin = glm::vec3(node->origin) + glm::vec3(offsetDir) * glm::vec3(node->size);

		// Traverse from root to find node with matching origin and same size
		std::function<OctreeNode*(OctreeNode*)> find = [&](OctreeNode* current) -> OctreeNode* {
			if (!current) return nullptr;

			float sizeDiff = fabs(current->size - node->size);
			float dist = glm::distance(glm::vec3(current->origin), neighborOrigin);

			if (sizeDiff < EPSILON && dist < EPSILON && current->isLeaf)
				return current;

			if (!current->isLeaf) {
				for (auto& child : current->children) {
					if (child) {
						OctreeNode* found = find(child.get());
						if (found) return found;
					}
				}
			}

			return nullptr;
		};

		return find(root);
	}

	struct BorderUpdate {
		OctreeNode* neighborNode;
		glm::ivec3 neighborPos;
		float newValue;
	};

	std::vector<BorderUpdate> deferredBorderUpdates;

bool ModifyScalarFieldInRadius(ScalarField& field, glm::ivec3 pos, float radius, float strength, OctreeNode* node, uint64_t planetUuid) {
    bool modified = false;

    constexpr int RES = SCALAR_FIELD_SIZE; // 16 for a 17x17x17 field
    glm::vec3 corner = glm::vec3(node->origin) - glm::vec3(node->size * 0.5f);
    float voxelSize = node->size / static_cast<float>(RES);

    for (int x = -static_cast<int>(radius); x <= static_cast<int>(radius); ++x) {
        for (int y = -static_cast<int>(radius); y <= static_cast<int>(radius); ++y) {
            for (int z = -static_cast<int>(radius); z <= static_cast<int>(radius); ++z) {
                glm::ivec3 offset(x, y, z);
                glm::ivec3 targetPos = pos + offset;

                // Ensure targetPos is within bounds
                if (targetPos.x < 0 || targetPos.x > RES ||
                    targetPos.y < 0 || targetPos.y > RES ||
                    targetPos.z < 0 || targetPos.z > RES)
                    continue;

                float dist = glm::length(glm::vec3(offset));
                if (dist > radius)
                    continue;

                float falloff = 1.0f - (dist / radius);
                float delta = strength * falloff;

                // Modify the scalar field

                modified = true;
            	node->scalarField[targetPos] += delta;
                // Check if it's a border voxel
                bool onBorder =
                    targetPos.x == 0 || targetPos.x == RES ||
                    targetPos.y == 0 || targetPos.y == RES ||
                    targetPos.z == 0 || targetPos.z == RES;

            	bool outerBorder = targetPos.x == RES||
targetPos.y == RES||
targetPos.z == RES;
            	if (!outerBorder) {
            		//node->scalarField[targetPos] += delta;
            	}
//
                //if (!onBorder)
                //    continue;

                // Compute world position at voxel center
                glm::vec3 worldPos = corner + (glm::vec3(targetPos) + 0.5f) * voxelSize;

                // Retrieve neighbor nodes
                auto* root = &static_cast<BigPhysicalBodyScript*>(Scene::GetActiveScene()
                    ->GetComponent<CppScriptComponent>(planetUuid)->mScripts[0])->mRootOctreeNode;

            	auto applyToNeighbor = [&](glm::ivec3 offset, glm::ivec3 thisPos, glm::ivec3 neighborPos) {
            		OctreeNode* neighbor = GetNodeAtOffset(node, offset, root);
            		if (!neighbor) return;
            		//neighbor->scalarField[neighborPos] += delta;
            		//neighbor->scalarField[neighborPos] = node->scalarField[targetPos];
            		deferredBorderUpdates.push_back({ neighbor, neighborPos, node->scalarField[targetPos] });
            	};


            	if (targetPos.x == 0) {
            		//applyToNeighbor(glm::ivec3(-1, 0, 0), targetPos, glm::ivec3(RES, targetPos.y, targetPos.z));
            	}
            	if (targetPos.x == RES) {
            		applyToNeighbor(glm::ivec3(1, 0, 0), targetPos, glm::ivec3(0, targetPos.y, targetPos.z));
            	}

            	if (targetPos.y == 0) {
            		//applyToNeighbor(glm::ivec3(0, -1, 0), targetPos, glm::ivec3(targetPos.x, RES, targetPos.z));
            	}
            	if (targetPos.y == RES) {
            		applyToNeighbor(glm::ivec3(0, 1, 0), targetPos, glm::ivec3(targetPos.x, 0, targetPos.z));
            	}

            	if (targetPos.z == 0) {
            		//applyToNeighbor(glm::ivec3(0, 0, -1), targetPos, glm::ivec3(targetPos.x, targetPos.y, RES));
            	}
            	if (targetPos.z == RES) {
            		applyToNeighbor(glm::ivec3(0, 0, 1), targetPos, glm::ivec3(targetPos.x, targetPos.y, 0));
            	}
                //std::vector<OctreeNode*> neighbors = GetNeighborNodes(node, root);
                //for (OctreeNode* neighbor : neighbors) {
                //    if (!neighbor) continue;
//
                //    glm::vec3 neighborCorner = neighbor->origin - glm::vec3(neighbor->size * 0.5f);
                //    float neighborVoxelSize = neighbor->size / static_cast<float>(RES);
//
                //    // Convert world position to neighbor-local voxel index
                //    glm::vec3 neighborLocal = (worldPos - neighborCorner) / neighborVoxelSize;
                //    glm::ivec3 neighborPos = glm::floor(neighborLocal);
//
                //    // Clamp neighborPos to valid range
                //    neighborPos = glm::clamp(neighborPos, glm::ivec3(0), glm::ivec3(RES));
//
                //    // Synchronize scalar field value
                //    neighbor->scalarField[neighborPos] += delta * 0.25f;//node->scalarField[targetPos];
                //}
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

		glm::vec3 min = glm::vec3(node->origin) - glm::vec3(node->size * 0.5f);
		glm::vec3 max = glm::vec3(node->origin) + glm::vec3(node->size * 0.5f);

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
				glm::vec3 corner = glm::vec3(hitNode->origin) - glm::vec3(hitNode->size * 0.5f);
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
				//affectedNodes.push_back(hitNode);
				CollectIntersectingLeafNodes(
					&static_cast<BigPhysicalBodyScript*>(scene->GetComponent<CppScriptComponent>(planet->uuid)->mScripts[0])->mRootOctreeNode,
					worldPosition,
					mDigRadius,
					affectedNodes);

				std::vector<OctreeNode*> nodesToUpdate;
				for (OctreeNode* node : affectedNodes) {
					// Transform worldPos into this node scalar field space
					glm::vec3 corner = glm::vec3(node->origin) - glm::vec3(node->size * 0.5f);
					float voxelSize = node->size / float(SCALAR_FIELD_SIZE);
					glm::vec3 local = (worldPosition - corner) / voxelSize;
					glm::ivec3 scalarFieldPos = glm::clamp(glm::ivec3(local), glm::ivec3(0), glm::ivec3(SCALAR_FIELD_SIZE));

					bool modified = ModifyScalarFieldInRadius(node->scalarField, scalarFieldPos, mDigRadius / (node->size / SCALAR_FIELD_SIZE), strength, node, planet->uuid);

					if (modified) {
						nodesToUpdate.push_back(node);
					}
				}

				for (const auto& update : deferredBorderUpdates) {
					update.neighborNode->scalarField[update.neighborPos] = update.newValue;
				}
				deferredBorderUpdates.clear();

				for (OctreeNode* node : nodesToUpdate) {
					PlanetGenerator::UpdateNodeMesh(scene, node);
				}
			}
		}
	}
} // namespace Plaza
#include "PlanetGenerator.h"
#include "MarchingCubes.h"
#include "Physics/PhysicalBodyScript.h"
#include "Physics/BigPhysicalBodyScript.h"
#include <thread>

namespace Plaza {
	void PlanetGenerator::OnStart(Scene* scene) { ECS::RegisterComponents(); }

	void collectLeavesInSphere(OctreeNode& node, std::vector<OctreeNode*>& outLeaves, const glm::vec3& sphereCenter,
							   float sphereRadius) {
		if (node.isLeaf) {
			outLeaves.push_back(&node);
			return;
		}

		for (auto& child : node.children) {
			if (child)
				collectLeavesInSphere(*child.get(), outLeaves, sphereCenter, sphereRadius);
		}
	}

	ScalarField PlanetGenerator::GenerateSphere(const glm::vec3& origin, float size, float planetRadius,
												PerlinNoise& perlin) {
		const int resolution = 33;
		float voxelSize = size / (resolution - 1); // Adjusted to use resolution - 1
		ScalarField field = ScalarField();

		for (int z = 0; z < resolution; ++z) {
			for (int y = 0; y < resolution; ++y) {
				for (int x = 0; x < resolution; ++x) {
					glm::vec3 localPos = glm::vec3(x, y, z) * voxelSize - glm::vec3(size / 2.0f);
					glm::vec3 worldPos = origin + localPos;

					// Normalize position to unit sphere space
					glm::vec3 normPos = glm::normalize(worldPos);
					float baseDist = glm::length(worldPos);

					// Layered Perlin noise (like in your other functions)
					float noise = 0.0f;
					float amplitude = mNoiseAmplitude;
					float frequency = mNoiseFrequency;
					float lodScaleFactor = 1.0f; // Example scaling

					// Calculate scale factor based on chunk size
					float scaleFactor = size / mPlanetRadius; // Adjust this based on your system
					for (int i = 0; i < 6; ++i) {
						// Adjust frequency based on chunk size
						float currentFrequency = frequency * (0.5f + lodScaleFactor) * scaleFactor;
						float currentAmplitude = amplitude * (0.5f + (1.0f - lodScaleFactor)); // Adjust amplitude

						// Sample the noise at the appropriate scale
						float ridged = 1.0f - fabs(perlin.Noise(normPos * currentFrequency));
						noise += ridged * ridged * currentAmplitude;
						amplitude *= 0.5f;
						frequency *= 2.0f;
					}

					float distanceToSurface = baseDist - planetRadius + (noise * size / (mPlanetRadius * 2.0f));

					if (noise == 0.0f)
						PL_INFO("Noise is 0 lol");

					// Use noise to modulate the surface
					field[glm::ivec3(x, y, z)] = distanceToSurface;
				}
			}
		}
		return field;
	}

	void PlanetGenerator::OnUpdateEditorGUI(Scene* scene) {
		if (ImGui::Button("Generate Planet")) {
			mPlanetEntity = scene->NewEntity("Planet");

			// std::thread worker([&] {
			std::mutex mutex;
			OctreeNode root;
			root.origin = glm::vec3(0.0f);	  // Center of the planet
			root.size = mPlanetRadius * 2.0f; // Full diameter
			root.isLeaf = true;
			root.Update(glm::vec3(0.0f, mPlanetRadius * 1.0f, 0.0f), mPlanetRadius);

			std::vector<OctreeNode*> renderNodes;
			collectLeavesInSphere(root, renderNodes, root.origin, mPlanetRadius);
			PL_INFO("Nodes: {}", renderNodes.size());

			if (renderNodes.size() > 10000)
				return;

			std::for_each(std::execution::par, renderNodes.begin(), renderNodes.end(), [&](auto&& node) {
				if (node == nullptr)
					return;
				ScalarField field = ScalarField();
				if (node->scalarField.size() <= 0) {
					node->scalarField = GenerateSphere(node->origin, node->size + 1, mPlanetRadius, mPerlin);
				}

				Mesh* chunkMesh = GenerateMesh(node->scalarField, node->origin, node->size, mIsoLevel,
											   node->size / 32.0f, node->origin);

				if (chunkMesh) {
					std::lock_guard lock(mutex);
					Entity* chunkEntity = scene->NewEntity("PlanetChunk " + std::to_string(node->size), mPlanetEntity);
					node->entityUuid = chunkEntity->uuid;
					ECS::TransformSystem::SetLocalPosition(*scene->GetComponent<TransformComponent>(chunkEntity->uuid),
														   scene, node->origin - glm::vec3(node->size / 2.0f));
					ECS::TransformSystem::SetLocalScale(*scene->GetComponent<TransformComponent>(chunkEntity->uuid),
														scene, glm::vec3(1.0f));

					MeshRenderer* meshRenderer = scene->NewComponent<MeshRenderer>(chunkEntity->uuid);
					meshRenderer->ChangeMesh(chunkMesh);
					meshRenderer->AddMaterial(AssetsManager::GetDefaultMaterial());

					Collider* collider = scene->NewComponent<Collider>(chunkEntity->uuid);
					collider->AddMeshShape(chunkMesh);
					ECS::ColliderSystem::InitCollider(
						collider, scene->GetComponent<TransformComponent>(chunkEntity->uuid), nullptr);
				}
			});
			scene->NewComponent<CppScriptComponent>(mPlanetEntity->uuid);
			auto* bigBody = static_cast<BigPhysicalBodyScript*>(
				scene->GetComponent<CppScriptComponent>(mPlanetEntity->uuid)
					->AddScriptNewInstance(scene,
										   AssetsManager::GetScriptByName("BigPhysicalBodyScript.h")->mAssetUuid));

			// Get gravity values for the body
			double g = 9.81; // Desired gravity at surface
			double planetMass = g * (mPlanetRadius * mPlanetRadius) / 6.67430e-11;
			double planetVolume = (4.0 / 3.0) * glm::pi<double>() * pow(mPlanetRadius, 3);
			double planetDensity = planetMass / planetVolume;

			bigBody->mRadius = mPlanetRadius;
			bigBody->mWeight = planetMass;
			bigBody->mPosition = glm::vec3(0.0f);
			bigBody->mRotation = glm::quat();
			bigBody->mRootOctreeNode = std::move(root);

			PhysicalBodyScript::AddValueToPhysicalBodies(mPlanetEntity->uuid, bigBody);
			//});
		}
		ImGui::DragInt("Grid Size", &mGridSize);
		ImGui::DragFloat("Radius", &mPlanetRadius);
		ImGui::DragFloat("Iso Level", &mIsoLevel);
		ImGui::DragFloat("Noise Frequency", &mNoiseFrequency);
		ImGui::DragFloat("Noise Amplitude", &mNoiseAmplitude);
		ImGui::DragFloat("Max Height", &mMaxHeight);
		ImGui::DragFloat("Min Height", &mMinHeight);
		ImGui::DragInt("Scale", &mScale);
	}

	Mesh* PlanetGenerator::GenerateMesh(const std::unordered_map<glm::ivec3, float, IVec3Comparator>& grid,
										const glm::vec3& offset, const int gridSize, const float isoLevel,
										const float resolution, const glm::ivec3 chunkOffset) {
		std::vector<glm::vec3> vertices, normals, tangents;
		std::vector<glm::vec2> uvs;
		std::vector<unsigned int> indices, materials{0};

		if (grid.size() <= 0)
			return nullptr;

		int size = 33;
		for (int x = 0; x < size - 1; ++x) {
			for (int y = 0; y < size - 1; ++y) {
				for (int z = 0; z < size - 1; ++z) {
					float cube[8] = {grid.at({x, y, z}),
									 grid.at({x + 1, y, z}),
									 grid.at({x + 1, y, z + 1}),
									 grid.at({x, y, z + 1}),
									 grid.at({x, y + 1, z}),
									 grid.at({x + 1, y + 1, z}),
									 grid.at({x + 1, y + 1, z + 1}),
									 grid.at({x, y + 1, z + 1})};
					glm::vec3 pos(x, y, z);
					MarchingCubes::March(resolution, cube, isoLevel, vertices, indices,
						normals, uvs, pos, gridSize,false);
				}
			}
		}

		if (vertices.size() <= 3)
			return nullptr;
		Mesh* mesh = Application::Get()->mRenderer->CreateNewMesh(vertices, normals, uvs, tangents, indices, materials,
																  false, {}, {});
		AssetsManager::AddMesh(mesh);
		return mesh;
	}

	void PlanetGenerator::UpdateNodeMesh(Scene* scene, OctreeNode* node) {
		Mesh* chunkMesh = GenerateMesh(node->scalarField, node->origin, node->size, 0.25f,
							   node->size / 32.0f, node->origin);

		if (chunkMesh) {
			MeshRenderer* meshRenderer = scene->GetComponent<MeshRenderer>(node->entityUuid);
			meshRenderer->ChangeMesh(chunkMesh);

			Collider* collider = scene->GetComponent<Collider>(node->entityUuid);
			collider->mShapes.clear();
			collider->AddMeshShape(chunkMesh);
			ECS::ColliderSystem::InitCollider(
				collider, scene->GetComponent<TransformComponent>(node->entityUuid), nullptr);
		}
	}

	const glm::vec3 offsets[8] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
								  {-1, -1, 1},	{1, -1, 1},	 {-1, 1, 1},  {1, 1, 1}};

	void OctreeNode::SubDivide() {
		float childSize = size / 2.0f;

		for (int i = 0; i < 8; ++i) {
			glm::vec3 offset = offsets[i] * (childSize / 2.0f);
			glm::vec3 childOrigin = origin + offset;

			children[i] = std::make_unique<OctreeNode>();
			children[i]->origin = childOrigin;
			children[i]->size = childSize;
			children[i]->isLeaf = true;
		}

		isLeaf = false;
	}

	void OctreeNode::Update(const glm::vec3& playerPosition, const float planetRadius) {
		// if (IsCompletelyOutsideSphere(glm::vec3(0.0f), planetRadius)) {
		//	// Don't subdivide or update
		//	PL_INFO("Outside sphere");
		//	return;
		// }
		if (glm::length(origin) - size * 0.5f > planetRadius)
			return;

		glm::vec3 halfSize = glm::vec3(size / 2.0f);
		glm::vec3 minBounds = origin - halfSize;
		glm::vec3 maxBounds = origin + halfSize;

		bool inside = playerPosition.x >= minBounds.x && playerPosition.x <= maxBounds.x &&
					  playerPosition.y >= minBounds.y && playerPosition.y <= maxBounds.y &&
					  playerPosition.z >= minBounds.z && playerPosition.z <= maxBounds.z;

		if (size == planetRadius)
			inside = true;

		float minSize = 32.0f;
		if (inside && size > minSize) {
			if (isLeaf)
				SubDivide();
			for (auto& child : children) {
				if (child)
					child->Update(playerPosition, planetRadius);
			}
			isLeaf = false;
		}
		else if (!inside || size <= minSize) {
			// Too small for current distance — collapse node
			for (auto& child : children) {
				child.reset();
			}
			isLeaf = true;
		}
	}

	bool OctreeNode::IsCompletelyOutsideSphere(const glm::vec3& center, float radius) {
		glm::vec3 halfSize = glm::vec3(size / 2.0f);
		glm::vec3 minBounds = origin - halfSize;
		glm::vec3 maxBounds = origin + halfSize;

		// Clamp center to AABB
		glm::vec3 closestPoint = glm::clamp(center, minBounds, maxBounds);
		float distance = glm::distance(center, closestPoint);

		return distance > radius;
	}
} // namespace Plaza

#include "PlanetGenerator.h"

#include "Constants.h"
#include "MarchingCubes.h"
#include "Physics/PhysicalBodyScript.h"
#include "Physics/BigPhysicalBodyScript.h"
#include <thread>

namespace Plaza {
	void PlanetGenerator::OnStart(Scene* scene) {
		ECS::RegisterComponents();
		mMaterialsVector =
			scene->GetComponent<MeshRenderer>(CppHelper::FindEntity("MaterialsEntity")->uuid)->GetMaterials();
	}

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
		const int resolution = SCALAR_FIELD_SIZE_WITH_PADDING;
		float voxelSize = size / (SCALAR_FIELD_SIZE);
		ScalarField field = ScalarField();

		float lodLevel = glm::log2(planetRadius / size);
		float frequencyLodMultiplier = std::pow(2.0f, lodLevel);
		float amplitudeLodMultiplier = std::pow(0.5f, lodLevel);

		for (int z = 0; z < resolution; ++z) {
			for (int y = 0; y < resolution; ++y) {
				for (int x = 0; x < resolution; ++x) {
					glm::ivec3 localPos = glm::ivec3(x, y, z) * glm::ivec3(voxelSize) - glm::ivec3(size * 0.5f);
					glm::vec3 worldPos = origin + glm::vec3(localPos);

					glm::vec3 normPos = glm::normalize(worldPos);
					float baseDist = glm::length(worldPos);

					localPos = glm::ivec3(x, y, z);
					glm::vec3 worldOffset = glm::vec3(localPos) * voxelSize;
					worldPos = origin + worldOffset;

					float noise = 0.0f;
					float amplitude = mNoiseAmplitude * amplitudeLodMultiplier;
					float frequency = mNoiseFrequency * frequencyLodMultiplier;
					float lodScaleFactor = glm::round(SCALAR_FIELD_SIZE / (size));

					float scaleFactor = size / mPlanetRadius;

					for (int i = 0; i < 6; ++i) {
						float ridged = 1.0f - fabs(perlin.Noise(worldPos * frequency));
						noise += ridged * ridged * amplitude;

						amplitude *= 0.5f;
						frequency *= 2.0f;
					}

					float distanceToSurface = (baseDist - planetRadius) + (noise * mMaxHeight);

					if (noise == 0.0f)
						PL_INFO("Noise is 0 lol");

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
			root.Update(glm::vec3(0.0f, mPlanetRadius * 1.0f, 0.0f), mPlanetRadius, mGenerateSmallestLOD);

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
					node->scalarField = GenerateSphere(node->origin, node->size, mPlanetRadius, mPerlin);
				}

				Mesh* chunkMesh = GenerateMesh(node->scalarField, node->origin, node->size, mIsoLevel,
											   node->size / SCALAR_FIELD_SIZE, node->origin);

				if (chunkMesh) {
					std::lock_guard lock(mutex);
					std::string originString = std::to_string(node->origin.x) + std::to_string(node->origin.y) +
											   std::to_string(node->origin.z);
					Entity* chunkEntity =
						scene->NewEntity("PlanetChunk " + std::to_string(node->size) + originString, mPlanetEntity);
					node->entityUuid = chunkEntity->uuid;
					ECS::TransformSystem::SetLocalPosition(*scene->GetComponent<TransformComponent>(chunkEntity->uuid),
														   scene,
														   glm::vec3(node->origin) - glm::vec3(node->size * 0.5f));
					ECS::TransformSystem::SetLocalScale(*scene->GetComponent<TransformComponent>(chunkEntity->uuid),
														scene, glm::vec3(1.0f));

					MeshRenderer* meshRenderer = scene->NewComponent<MeshRenderer>(chunkEntity->uuid);
					meshRenderer->ChangeMesh(chunkMesh);
					for (Material* material : mMaterialsVector) {
						meshRenderer->AddMaterial(material);
					}

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
		ImGui::DragInt("Radius", &mPlanetRadius);
		ImGui::DragFloat("Iso Level", &mIsoLevel);
		ImGui::DragFloat("Noise Frequency", &mNoiseFrequency, 0.01f, 0.0f, 1000.0f, "%.6f");
		ImGui::DragFloat("Noise Amplitude", &mNoiseAmplitude);
		ImGui::DragFloat("Max Height", &mMaxHeight);
		ImGui::DragFloat("Min Height", &mMinHeight);
		ImGui::DragInt("Scale", &mScale);
		ImGui::Checkbox("Generate Smallest LOD", &mGenerateSmallestLOD);
	}

	Mesh* PlanetGenerator::GenerateMesh(const std::unordered_map<glm::ivec3, float, IVec3Comparator>& grid,
										const glm::vec3& offset, const int gridSize, const float isoLevel,
										const float resolution, const glm::ivec3 chunkOffset) {
		std::vector<glm::vec3> vertices, normals, tangents;
		std::vector<glm::vec2> uvs;
		std::vector<unsigned int> indices, materials{0};

		if (grid.size() <= 0)
			return nullptr;

		size_t estimatedMaxTriangles = SCALAR_FIELD_SIZE * SCALAR_FIELD_SIZE * SCALAR_FIELD_SIZE;
		vertices.reserve(estimatedMaxTriangles * 3);
		normals.reserve(estimatedMaxTriangles * 3);
		uvs.reserve(estimatedMaxTriangles * 3);
		indices.reserve(estimatedMaxTriangles * 3);
		materials.reserve(estimatedMaxTriangles * 3);

		for (int x = 0; x < SCALAR_FIELD_SIZE; ++x) {
			for (int y = 0; y < SCALAR_FIELD_SIZE; ++y) {
				for (int z = 0; z < SCALAR_FIELD_SIZE; ++z) {
					float cube[8] = {grid.at({x, y, z}),
									 grid.at({x + 1, y, z}),
									 grid.at({x + 1, y, z + 1}),
									 grid.at({x, y, z + 1}),
									 grid.at({x, y + 1, z}),
									 grid.at({x + 1, y + 1, z}),
									 grid.at({x + 1, y + 1, z + 1}),
									 grid.at({x, y + 1, z + 1})};
					glm::vec3 pos(x, y, z);
					MarchingCubes::March(resolution, cube, isoLevel, vertices, indices, normals, uvs, materials, pos, gridSize,
										 false);
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
		Mesh* chunkMesh = GenerateMesh(node->scalarField, node->origin, node->size, 0.0f,
									   node->size / SCALAR_FIELD_SIZE, node->origin);

		if (chunkMesh) {
			MeshRenderer* meshRenderer = scene->GetComponent<MeshRenderer>(node->entityUuid);
			if (!meshRenderer)
				return;
			meshRenderer->ChangeMesh(chunkMesh);

			Collider* collider = scene->GetComponent<Collider>(node->entityUuid);
			collider->mShapes.clear();
			collider->AddMeshShape(chunkMesh);
			ECS::ColliderSystem::InitCollider(collider, scene->GetComponent<TransformComponent>(node->entityUuid),
											  nullptr);
		}
	}
} // namespace Plaza

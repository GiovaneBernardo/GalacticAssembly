#include "PlanetGenerator.h"

#include "Constants.h"
#include "MarchingCubes.h"
#include "Physics/PhysicalBodyScript.h"
#include "Physics/BigPhysicalBodyScript.h"
#include <thread>

namespace Plaza {

	// Transvoxel face directions
	enum class TransvoxelFace : int {
		NEGATIVE_X = 0,
		POSITIVE_X = 1,
		NEGATIVE_Y = 2,
		POSITIVE_Y = 3,
		NEGATIVE_Z = 4,
		POSITIVE_Z = 5
	};

	// Enhanced OctreeNode with neighbor tracking
	struct EnhancedOctreeNode : public OctreeNode {
		std::array<EnhancedOctreeNode*, 6> neighbors{nullptr}; // 6 face neighbors
		std::array<bool, 6> needsTransition{false}; // Which faces need transition cells
		int lodLevel = 0;

		void SubDivide() override {
			float childSize = size / 2.0f;
			const glm::vec3 offsets[8] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
										  {-1, -1, 1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};
			for (int i = 0; i < 8; ++i) {
				glm::vec3 offset = offsets[i];
				glm::vec3 childOrigin = glm::vec3(origin) + offset * (childSize / 2.0f);
				float voxelSize = childSize / SCALAR_FIELD_SIZE;
				childOrigin = glm::round(childOrigin / voxelSize) * voxelSize;

				auto child = std::make_unique<EnhancedOctreeNode>();
				child->origin = childOrigin;
				child->size = childSize;
				child->isLeaf = true;
				children[i] = std::move(child);
			}
			isLeaf = false;
		}

		bool HasLowerLODNeighbor(TransvoxelFace face) const {
			EnhancedOctreeNode* neighbor = neighbors[static_cast<int>(face)];
			return neighbor && (neighbor->lodLevel < this->lodLevel);
		}

		void UpdateNeighborInfo() {
			for (int i = 0; i < 6; ++i) {
				needsTransition[i] = HasLowerLODNeighbor(static_cast<TransvoxelFace>(i));
			}
		}
	};

	void PlanetGenerator::OnStart(Scene* scene) {
		ECS::RegisterComponents();
		mMaterialsVector =
			scene->GetComponent<MeshRenderer>(CppHelper::FindEntity("MaterialsEntity")->uuid)->GetMaterials();

		GeneratePlanet(scene);
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

	// Build neighbor relationships in octree
	void PlanetGenerator::BuildNeighborRelationships(EnhancedOctreeNode& root) {
		std::vector<EnhancedOctreeNode*> allNodes;
		CollectAllNodes(root, allNodes);

		for (EnhancedOctreeNode* node : allNodes) {
			if (!node->isLeaf) continue;

			// Find neighbors for each face
			for (int face = 0; face < 6; ++face) {
				node->neighbors[face] = FindNeighborAtFace(node, static_cast<TransvoxelFace>(face), allNodes);
			}

			// Calculate LOD level based on size
			node->lodLevel = static_cast<int>(glm::log2(static_cast<float>(root.size) / static_cast<float>(node->size)));
			node->UpdateNeighborInfo();
		}
	}

	void PlanetGenerator::CollectAllNodes(EnhancedOctreeNode& node, std::vector<EnhancedOctreeNode*>& outNodes) {
		outNodes.push_back(&node);

		if (!node.isLeaf) {
			for (auto& child : node.children) {
				if (child) {
					CollectAllNodes(*static_cast<EnhancedOctreeNode*>(child.get()), outNodes);
				}
			}
		}
	}

	EnhancedOctreeNode* PlanetGenerator::FindNeighborAtFace(EnhancedOctreeNode* node, TransvoxelFace face,
															const std::vector<EnhancedOctreeNode*>& allNodes) {
		glm::vec3 faceOffset = GetFaceOffset(face) * (node->size * 0.5f + 0.1f); // Small epsilon
		glm::vec3 searchPos = glm::vec3(node->origin) + faceOffset;

		EnhancedOctreeNode* bestCandidate = nullptr;
		float bestDistance = std::numeric_limits<float>::max();

		for (EnhancedOctreeNode* candidate : allNodes) {
			if (candidate == node || !candidate->isLeaf) continue;

			// Check if search position is within candidate bounds
			glm::vec3 halfSize = glm::vec3(candidate->size * 0.5f);
			glm::vec3 minBounds = glm::vec3(candidate->origin) - halfSize;
			glm::vec3 maxBounds = glm::vec3(candidate->origin) + halfSize;

			if (searchPos.x >= minBounds.x && searchPos.x <= maxBounds.x &&
				searchPos.y >= minBounds.y && searchPos.y <= maxBounds.y &&
				searchPos.z >= minBounds.z && searchPos.z <= maxBounds.z) {

				float distance = glm::distance(searchPos, glm::vec3(candidate->origin));
				if (distance < bestDistance) {
					bestDistance = distance;
					bestCandidate = candidate;
				}
			}
		}

		return bestCandidate;
	}

	glm::vec3 PlanetGenerator::GetFaceOffset(TransvoxelFace face) {
		switch (face) {
			case TransvoxelFace::NEGATIVE_X: return glm::vec3(-1, 0, 0);
			case TransvoxelFace::POSITIVE_X: return glm::vec3(1, 0, 0);
			case TransvoxelFace::NEGATIVE_Y: return glm::vec3(0, -1, 0);
			case TransvoxelFace::POSITIVE_Y: return glm::vec3(0, 1, 0);
			case TransvoxelFace::NEGATIVE_Z: return glm::vec3(0, 0, -1);
			case TransvoxelFace::POSITIVE_Z: return glm::vec3(0, 0, 1);
			default: return glm::vec3(0);
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

					for (int i = 0; i < 6; ++i) {
						float ridged = 1.0f - fabs(perlin.Noise(worldPos * frequency));
						noise += ridged * ridged * amplitude;

						amplitude *= 0.5f;
						frequency *= 2.0f;
					}

					float distanceToSurface = (baseDist - planetRadius) + (noise * mMaxHeight);

					field[glm::ivec3(x, y, z)] = distanceToSurface;
				}
			}
		}
		return field;
	}

	void PlanetGenerator::OnUpdateEditorGUI(Scene* scene) {
		if (ImGui::Button("Generate Planet")) {
			GeneratePlanet(scene);
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

	Entity* PlanetGenerator::GeneratePlanet(Scene* scene) {
		mPlanetEntity = scene->NewEntity("Planet");

		// Use enhanced octree node
		auto enhancedRoot = new EnhancedOctreeNode();
		enhancedRoot->origin = glm::vec3(0.0f);
		enhancedRoot->size = mPlanetRadius * 2.0f;
		enhancedRoot->isLeaf = true;
		enhancedRoot->Update(glm::vec3(0.0f, mPlanetRadius * 1.0f, 0.0f), mPlanetRadius, mGenerateSmallestLOD);

		// Build neighbor relationships
		BuildNeighborRelationships(*enhancedRoot);

		std::vector<OctreeNode*> renderNodes;
		collectLeavesInSphere(*enhancedRoot, renderNodes, enhancedRoot->origin, mPlanetRadius);
		PL_INFO("Nodes: {}", renderNodes.size());

		if (renderNodes.size() > 10000)
			return nullptr;

		// Phase 1: Generate scalar fields in parallel (CPU-heavy, thread-safe)
		std::for_each(std::execution::par, renderNodes.begin(), renderNodes.end(), [&](auto&& node) {
			if (node == nullptr)
				return;

			EnhancedOctreeNode* enhancedNode = static_cast<EnhancedOctreeNode*>(node);

			if (enhancedNode->scalarField.size() <= 0) {
				enhancedNode->scalarField =
					GenerateSphere(enhancedNode->origin, enhancedNode->size, mPlanetRadius, mPerlin);
			}
		});

		// Phase 2: Create meshes, entities, and physics sequentially
		// (CreateNewMesh, AssetsManager, PhysX cooking are not thread-safe)
		for (auto& node : renderNodes) {
			if (node == nullptr)
				continue;

			EnhancedOctreeNode* enhancedNode = static_cast<EnhancedOctreeNode*>(node);
			Mesh* chunkMesh = GenerateTransvoxelMesh(enhancedNode, mIsoLevel);

			if (chunkMesh) {
				std::string originString = std::to_string(enhancedNode->origin.x) +
										   std::to_string(enhancedNode->origin.y) +
										   std::to_string(enhancedNode->origin.z);
				Entity* chunkEntity =
					scene->NewEntity("PlanetChunk " + std::to_string(enhancedNode->size) + originString, mPlanetEntity);
				enhancedNode->entityUuid = chunkEntity->uuid;
				ECS::TransformSystem::SetLocalPosition(
					*scene->GetComponent<TransformComponent>(chunkEntity->uuid), scene,
					glm::vec3(enhancedNode->origin) - glm::vec3(enhancedNode->size * 0.5f));
				ECS::TransformSystem::SetLocalScale(*scene->GetComponent<TransformComponent>(chunkEntity->uuid), scene,
													glm::vec3(1.0f));

				MeshRenderer* meshRenderer = scene->NewComponent<MeshRenderer>(chunkEntity->uuid);
				meshRenderer->ChangeMesh(chunkMesh);
				for (Material* material : mMaterialsVector) {
					meshRenderer->AddMaterial(material);
				}

				Collider* collider = scene->NewComponent<Collider>(chunkEntity->uuid);
				collider->AddMeshShape(chunkMesh);
				ECS::ColliderSystem::InitCollider(collider, scene->GetComponent<TransformComponent>(chunkEntity->uuid),
												  nullptr);
			}
		}

		scene->NewComponent<CppScriptComponent>(mPlanetEntity->uuid);
		auto* bigBody = static_cast<BigPhysicalBodyScript*>(
			scene->GetComponent<CppScriptComponent>(mPlanetEntity->uuid)
				->AddScriptNewInstance(scene, AssetsManager::GetScriptByName("BigPhysicalBodyScript.h")->mAssetUuid));

		double g = 9.81;
		double planetMass = g * (mPlanetRadius * mPlanetRadius) / 6.67430e-11;
		double planetVolume = (4.0 / 3.0) * glm::pi<double>() * pow(mPlanetRadius, 3);

		bigBody->mRadius = mPlanetRadius;
		bigBody->mWeight = planetMass;
		bigBody->mPosition = glm::vec3(0.0f);
		bigBody->mRotation = glm::quat();
		// Move the octree data to BigPhysicalBodyScript
		//bigBody->mRootOctreeNode = std::move(static_cast<OctreeNode&>(*enhancedRoot));
		bigBody->mRootOctreeNode = std::move(*enhancedRoot);
		return mPlanetEntity;
	}

	// Enhanced mesh generation with transvoxel support
	Mesh* PlanetGenerator::GenerateTransvoxelMesh(EnhancedOctreeNode* node, float isoLevel) {
		std::vector<glm::vec3> vertices, normals, tangents;
		std::vector<glm::vec2> uvs;
		std::vector<unsigned int> indices;
		std::vector<unsigned int> materials;

		const ScalarField& grid = node->scalarField;
		if (grid.size() <= 0)
			return nullptr;

		size_t estimatedMaxTriangles = SCALAR_FIELD_SIZE * SCALAR_FIELD_SIZE * SCALAR_FIELD_SIZE;
		vertices.reserve(estimatedMaxTriangles * 3);
		normals.reserve(estimatedMaxTriangles * 3);
		materials.reserve(estimatedMaxTriangles * 3);
		uvs.reserve(estimatedMaxTriangles * 3);
		indices.reserve(estimatedMaxTriangles * 3);

		float resolution = node->size / SCALAR_FIELD_SIZE;

		// 1. Generate regular marching cubes for interior
		GenerateInteriorMesh(grid, node, isoLevel, resolution, vertices, indices, normals, uvs, materials);

		// 2. Generate transition cells for faces that need them
		for (int face = 0; face < 6; ++face) {
			if (node->needsTransition[face]) {
				TransvoxelFace faceEnum = static_cast<TransvoxelFace>(face);
				EnhancedOctreeNode* neighbor = node->neighbors[face];

				if (neighbor && neighbor->scalarField.size() > 0) {
					GenerateTransitionCells(grid, neighbor->scalarField, faceEnum, node, neighbor,
										  isoLevel, resolution, vertices, indices, normals, uvs, materials);
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

	void PlanetGenerator::GenerateInteriorMesh(const ScalarField& grid, EnhancedOctreeNode* node, float isoLevel,
											  float resolution, std::vector<glm::vec3>& vertices,
											  std::vector<unsigned int>& indices, std::vector<glm::vec3>& normals,
											  std::vector<glm::vec2>& uvs, std::vector<unsigned int>& materials) {

		// Use existing marching cubes but avoid boundary cells that need transitions
		for (int x = 0; x < SCALAR_FIELD_SIZE; ++x) {
			for (int y = 0; y < SCALAR_FIELD_SIZE; ++y) {
				for (int z = 0; z < SCALAR_FIELD_SIZE; ++z) {

					// Skip boundary cells that will be handled by transition cells
					bool onBoundary = false;
					if ((x == 0 && node->needsTransition[static_cast<int>(TransvoxelFace::NEGATIVE_X)]) ||
						(x == SCALAR_FIELD_SIZE - 1 && node->needsTransition[static_cast<int>(TransvoxelFace::POSITIVE_X)]) ||
						(y == 0 && node->needsTransition[static_cast<int>(TransvoxelFace::NEGATIVE_Y)]) ||
						(y == SCALAR_FIELD_SIZE - 1 && node->needsTransition[static_cast<int>(TransvoxelFace::POSITIVE_Y)]) ||
						(z == 0 && node->needsTransition[static_cast<int>(TransvoxelFace::NEGATIVE_Z)]) ||
						(z == SCALAR_FIELD_SIZE - 1 && node->needsTransition[static_cast<int>(TransvoxelFace::POSITIVE_Z)])) {
						onBoundary = true;
					}

					if (onBoundary) continue;

					float cube[8] = {
						grid.at({x, y, z}),
						grid.at({x + 1, y, z}),
						grid.at({x + 1, y, z + 1}),
						grid.at({x, y, z + 1}),
						grid.at({x, y + 1, z}),
						grid.at({x + 1, y + 1, z}),
						grid.at({x + 1, y + 1, z + 1}),
						grid.at({x, y + 1, z + 1})
					};

					glm::vec3 pos(x, y, z);
					MarchingCubes::March(resolution, cube, isoLevel, vertices, indices, normals, uvs, materials, pos, node->size, false);
				}
			}
		}
	}

	void PlanetGenerator::GenerateTransitionCells(const ScalarField& highResField, const ScalarField& lowResField,
												 TransvoxelFace face, EnhancedOctreeNode* highResNode,
												 EnhancedOctreeNode* lowResNode, float isoLevel, float resolution,
												 std::vector<glm::vec3>& vertices, std::vector<unsigned int>& indices,
												 std::vector<glm::vec3>& normals, std::vector<glm::vec2>& uvs,
												 std::vector<unsigned int>& materials) {

		// Simplified transition cell generation
		// This is a basic implementation - full Transvoxel would use lookup tables

		int faceSize = SCALAR_FIELD_SIZE;

		// Get face coordinates based on face direction
		for (int v = 0; v < faceSize - 1; ++v) {
			for (int u = 0; u < faceSize - 1; ++u) {
				// Sample transition cell (simplified - would use 13 samples in full Transvoxel)
				float samples[4];
				glm::ivec3 coords[4];

				GetFaceCoordinates(face, u, v, coords);

				// Sample high-res field
				for (int i = 0; i < 4; ++i) {
					auto it = highResField.find(coords[i]);
					samples[i] = (it != highResField.end()) ? it->second : 0.0f;
				}

				// Generate transition geometry (simplified approach)
				GenerateTransitionQuad(samples, coords, face, highResNode, isoLevel, resolution,
									 vertices, indices, normals, uvs, materials);
			}
		}
	}

	void PlanetGenerator::GetFaceCoordinates(TransvoxelFace face, int u, int v, glm::ivec3 coords[4]) {
		switch (face) {
			case TransvoxelFace::NEGATIVE_X:
				coords[0] = {0, u, v};
				coords[1] = {0, u+1, v};
				coords[2] = {0, u+1, v+1};
				coords[3] = {0, u, v+1};
				break;
			case TransvoxelFace::POSITIVE_X:
				coords[0] = {SCALAR_FIELD_SIZE-1, u, v};
				coords[1] = {SCALAR_FIELD_SIZE-1, u+1, v};
				coords[2] = {SCALAR_FIELD_SIZE-1, u+1, v+1};
				coords[3] = {SCALAR_FIELD_SIZE-1, u, v+1};
				break;
			case TransvoxelFace::NEGATIVE_Y:
				coords[0] = {u, 0, v};
				coords[1] = {u+1, 0, v};
				coords[2] = {u+1, 0, v+1};
				coords[3] = {u, 0, v+1};
				break;
			case TransvoxelFace::POSITIVE_Y:
				coords[0] = {u, SCALAR_FIELD_SIZE-1, v};
				coords[1] = {u+1, SCALAR_FIELD_SIZE-1, v};
				coords[2] = {u+1, SCALAR_FIELD_SIZE-1, v+1};
				coords[3] = {u, SCALAR_FIELD_SIZE-1, v+1};
				break;
			case TransvoxelFace::NEGATIVE_Z:
				coords[0] = {u, v, 0};
				coords[1] = {u+1, v, 0};
				coords[2] = {u+1, v+1, 0};
				coords[3] = {u, v+1, 0};
				break;
			case TransvoxelFace::POSITIVE_Z:
				coords[0] = {u, v, SCALAR_FIELD_SIZE-1};
				coords[1] = {u+1, v, SCALAR_FIELD_SIZE-1};
				coords[2] = {u+1, v+1, SCALAR_FIELD_SIZE-1};
				coords[3] = {u, v+1, SCALAR_FIELD_SIZE-1};
				break;
		}
	}

	void PlanetGenerator::GenerateTransitionQuad(float samples[4], glm::ivec3 coords[4], TransvoxelFace face,
												EnhancedOctreeNode* node, float isoLevel, float resolution,
												std::vector<glm::vec3>& vertices, std::vector<unsigned int>& indices,
												std::vector<glm::vec3>& normals, std::vector<glm::vec2>& uvs,
												std::vector<unsigned int>& materials) {

		// Simplified transition quad generation
		// This creates basic geometry to prevent cracks - not full Transvoxel quality

		std::vector<glm::vec3> quadVertices;

		for (int i = 0; i < 4; ++i) {
			if ((samples[i] <= isoLevel && samples[(i+1)%4] > isoLevel) ||
				(samples[i] > isoLevel && samples[(i+1)%4] <= isoLevel)) {

				// Interpolate along edge
				float t = (isoLevel - samples[i]) / (samples[(i+1)%4] - samples[i]);
				glm::vec3 pos = glm::mix(glm::vec3(coords[i]), glm::vec3(coords[(i+1)%4]), t);

				// Transform to world space
				pos = pos * resolution + glm::vec3(node->origin) - glm::vec3(node->size * 0.5f);
				quadVertices.push_back(pos);
			}
		}

		// Create triangles from vertices (simplified)
		if (quadVertices.size() >= 3) {
			unsigned int baseIndex = vertices.size();

			for (const auto& vertex : quadVertices) {
				vertices.push_back(vertex);
				normals.push_back(glm::normalize(vertex)); // Simple normal calculation
				uvs.push_back({0.0f, 0.0f}); // Basic UV
			}

			// Create triangles (fan triangulation)
			for (size_t i = 1; i < quadVertices.size() - 1; ++i) {
				indices.push_back(baseIndex);
				indices.push_back(baseIndex + i);
				indices.push_back(baseIndex + i + 1);
			}

			materials.insert(materials.end(), quadVertices.size(), 0);
		}
	}

	// Keep original GenerateMesh for backward compatibility
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
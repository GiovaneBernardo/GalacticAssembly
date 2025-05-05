#include "PlanetGenerator.h"
#include "MarchingCubes.h"
#include "Physics/PhysicalBodyScript.h"
#include "Physics/BigPhysicalBodyScript.h"
#include <thread>

namespace Plaza {
    void PlanetGenerator::OnStart(Scene *scene) {
        ECS::RegisterComponents();
    }

    void collectLeavesInSphere(OctreeNode &node, std::vector<OctreeNode *> &outLeaves, const glm::vec3 &sphereCenter,
                               float sphereRadius) {
        if (node.isLeaf) {
            outLeaves.push_back(&node);
            return;
        }

        for (auto &child: node.children) {
            if (child) collectLeavesInSphere(*child.get(), outLeaves, sphereCenter, sphereRadius);
        }
    }

    ScalarField PlanetGenerator::GenerateSphere(const glm::vec3 &origin, float size, float planetRadius,
                                                PerlinNoise &perlin) {
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

                    for (int i = 0; i < 6; ++i) {
                        float currentFrequency = frequency * (0.5f + lodScaleFactor); // Adjust frequency based on LOD
                        float currentAmplitude = amplitude * (0.5f + (1.0f - lodScaleFactor)); // Adjust amplitude

                        float ridged = 1.0f - fabs(perlin.Noise(normPos * currentFrequency));
                        noise += ridged * ridged * currentAmplitude;
                        amplitude *= 0.5f;
                        frequency *= 2.0f;
                    }

                    float distanceToSurface = baseDist - planetRadius + (noise * size / (mPlanetRadius * 2.0f));

                    if (noise == 0.0f)
                        PL_INFO("Noise is 0 lol");

                    // Use noise to modulate the surface
                    //float distanceToSurface = baseDist - planetRadius + noise; //* (planetRadius * 0.01f);
                    field[glm::ivec3(x, y, z)] = distanceToSurface;
                }
            }
        }
        return field;
    }

    void PlanetGenerator::OnUpdateEditorGUI(Scene *scene) {
        if (ImGui::Button("Generate Planet")) {
            mPlanetEntity = scene->NewEntity("Planet");

           // std::thread worker([&] {
                OctreeNode root;
                root.origin = glm::vec3(0.0f); // Center of the planet
                root.size = mPlanetRadius * 2.0f; // Full diameter
                root.isLeaf = true;
                root.Update(glm::vec3(0.0f, mPlanetRadius * 1.0f, 0.0f), mPlanetRadius);

                std::vector<OctreeNode *> renderNodes;
                collectLeavesInSphere(root, renderNodes, root.origin, mPlanetRadius);
                PL_INFO("Nodes: {}", renderNodes.size());

                if (renderNodes.size() > 10000)
                    return;

                for (auto *node: renderNodes) {
                    if (node == nullptr)
                        continue;
                    ScalarField field = ScalarField();
                    if (node->scalarField.size() <= 0) {
                        node->scalarField = GenerateSphere(node->origin, node->size + 1, mPlanetRadius, mPerlin);
                    }

                    Mesh *chunkMesh = GenerateMesh(node->scalarField, node->origin, node->size, mIsoLevel,
                                                   node->size / 32.0f, node->origin);

                    if (chunkMesh) {
                        Entity *chunkEntity = scene->NewEntity("PlanetChunk " + std::to_string(node->size),
                                                               mPlanetEntity);
                        ECS::TransformSystem::SetLocalPosition(
                            *scene->GetComponent<TransformComponent>(chunkEntity->uuid), scene,
                            node->origin - glm::vec3(node->size / 2.0f));
                        ECS::TransformSystem::SetLocalScale(*scene->GetComponent<TransformComponent>(chunkEntity->uuid),
                                                            scene, glm::vec3(1.0f));

                        MeshRenderer *meshRenderer = scene->NewComponent<MeshRenderer>(chunkEntity->uuid);
                        meshRenderer->ChangeMesh(chunkMesh);
                        meshRenderer->AddMaterial(AssetsManager::GetDefaultMaterial());

                        Collider *collider = scene->NewComponent<Collider>(chunkEntity->uuid);
                        collider->AddMeshShape(chunkMesh);
                    }
                }
         //   });
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

    Mesh *PlanetGenerator::GenerateMesh(const std::unordered_map<glm::ivec3, float, IVec3Comparator> &grid,
                                        const glm::vec3 &offset,
                                        const int gridSize, const float isoLevel, const float resolution,
                                        const glm::ivec3 chunkOffset) {
        std::vector<glm::vec3> vertices, normals, tangents;
        std::vector<glm::vec2> uvs;
        std::vector<unsigned int> indices, materials{0};

        if (grid.size() <= 0)
            return nullptr;

        int size = 33;
        for (int x = 0; x < size - 1; ++x) {
            for (int y = 0; y < size - 1; ++y) {
                for (int z = 0; z < size - 1; ++z) {
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
                    MarchingCubes::March(resolution, cube, isoLevel, vertices, indices, normals, uvs, pos, gridSize,
                                         mNoiseFrequency, mMaxHeight, mMinHeight, false);
                }
            }
        }

        if (vertices.size() <= 3) return nullptr;
        Mesh *mesh = Application::Get()->mRenderer->CreateNewMesh(vertices, normals, uvs, tangents, indices, materials,
                                                                  false, {}, {});
        AssetsManager::AddMesh(mesh);
        return mesh;
    }

    void PlanetGenerator::GeneratePlanetChunks(Scene *scene, const int chunkSize, const int gridSize,
                                               const float isoLevel, const int scale) {
        Entity *planetEntity = scene->NewEntity("Planet");
        ECS::TransformSystem::SetLocalScale(*scene->GetComponent<TransformComponent>(planetEntity->uuid), scene,
                                            glm::vec3(scale));

        int numChunks = gridSize / chunkSize;
        glm::vec3 planetCenter = glm::vec3(gridSize / 2);

        std::unordered_map<glm::ivec3, BigPhysicalBodyScript::Chunk, IVec3Comparator> chunks;

        PerlinNoise noise;
        for (int x = 0; x < numChunks; ++x) {
            for (int y = 0; y < numChunks; ++y) {
                for (int z = 0; z < numChunks; ++z) {
                    glm::ivec3 chunkOffset = glm::ivec3(x, y, z) * chunkSize;
                    glm::vec3 chunkCenter = glm::vec3(chunkOffset) + glm::vec3(chunkSize / 2);

                    // Distance from camera or arbitrary point (improve later)
                    float distance = glm::distance(chunkCenter, glm::vec3(0.0f, gridSize, 0.0f));
                    int resolution = 32; //GetLODLevel(distance);

                    // Skip interior chunks for spherical shape
                    if (glm::length(chunkCenter - planetCenter) > gridSize)
                        continue;

                    PL_CORE_INFO("LOD: {}, Distance: {}, Chunk: ({}, {}, {})", resolution, distance, x, y, z);

                    std::unordered_map<glm::ivec3, float, IVec3Comparator> chunkGrid = GenerateChunkGrid(
                        noise, chunkOffset, resolution, gridSize, chunkSize, isoLevel, mNoiseFrequency, mMaxHeight);
                    chunks[glm::ivec3(x, y, z)] = BigPhysicalBodyScript::Chunk{chunkGrid};

                    Mesh *chunkMesh = GenerateMesh(chunkGrid, chunkOffset, resolution, isoLevel, resolution,
                                                   chunkOffset);

                    if (chunkMesh) {
                        Entity *chunkEntity = scene->NewEntity("PlanetChunk", planetEntity);
                        ECS::TransformSystem::SetLocalPosition(
                            *scene->GetComponent<TransformComponent>(chunkEntity->uuid), scene, chunkOffset);
                        //ECS::TransformSystem::SetLocalScale(*scene->GetComponent<TransformComponent>(chunkEntity->uuid), scene, glm::vec3((float)lodLevels[0].resolution / (float)resolution));

                        MeshRenderer *meshRenderer = scene->NewComponent<MeshRenderer>(chunkEntity->uuid);
                        meshRenderer->ChangeMesh(chunkMesh);
                        meshRenderer->AddMaterial(AssetsManager::GetDefaultMaterial());

                        Collider *collider = scene->NewComponent<Collider>(chunkEntity->uuid);
                        collider->AddMeshShape(chunkMesh);
                    }
                }
            }
        }

        scene->NewComponent<CppScriptComponent>(planetEntity->uuid);
        auto *bigBody = static_cast<BigPhysicalBodyScript *>(
            scene->GetComponent<CppScriptComponent>(planetEntity->uuid)->AddScriptNewInstance(
                scene, AssetsManager::GetScriptByName("BigPhysicalBodyScript.h")->mAssetUuid));

        bigBody->mRadius = gridSize * scale;
        bigBody->mWeight = gridSize * 100.0f;
        bigBody->mPosition = glm::vec3(0.0f);
        bigBody->mRotation = glm::quat();
        bigBody->mChunks = chunks;

        double G = 6.674e-11;
        double M = 5.972e24 * 0.01;
        PL_CORE_INFO("Gravity Estimate: {}", (G * M) / pow(gridSize * scale, 2));
        PhysicalBodyScript::mBigPhysicalBodies.emplace(planetEntity->uuid, bigBody);
    }

    std::unordered_map<glm::ivec3, float, IVec3Comparator> PlanetGenerator::GenerateChunkGrid(
        PerlinNoise &perlin, const glm::ivec3 &chunkOffset, const int resolution,
        const int planetRadius, const int chunkSize, const float isoLevel,
        const float noiseFrequency, const float noiseAmplitude) {
        std::unordered_map<glm::ivec3, float, IVec3Comparator> grid;
        float step = float(chunkSize) / resolution;

        float average = 0.0f;
        float highest = 0.0f;

        for (int x = 0; x <= resolution + 1; ++x)
            for (int y = 0; y <= resolution + 1; ++y)
                for (int z = 0; z <= resolution + 1; ++z) {
                    glm::vec3 localPos = glm::vec3(x, y, z) * step;
                    glm::vec3 global = glm::vec3(chunkOffset) + localPos;
                    glm::vec3 pos = (global - glm::vec3(planetRadius / 2.0f)) / float(planetRadius);

                    float base = glm::length(pos);
                    float noise = 0.0f, amp = noiseAmplitude, freq = noiseFrequency;

                    for (int i = 0; i < 6; ++i) {
                        float ridged = 1.0f - fabs(perlin.Noise(pos * freq));
                        noise += ridged * ridged * amp;
                        amp *= 0.5f;
                        freq *= 2.0f;
                    }

                    float surfaceDist = fabs(base - isoLevel);
                    float falloff = exp(-surfaceDist * 10.0f);
                    noise *= falloff;

                    //average += base - noise;
                    //highest = glm::max(highest, base - noise);
                    average += base;
                    highest = glm::max(highest, base);

                    //PL_INFO("Base: {}", base);
                    grid[{x, y, z}] = base; // - noise;
                }

        PL_INFO("Average: {}; Highest: {}", average / (resolution * resolution * resolution), highest);
        return grid;
    }

    std::unordered_map<glm::ivec3, float, IVec3Comparator> PlanetGenerator::GenerateGrid(
        const int gridSize, const float isoLevel, const float noiseFrequency, const float noiseAmplitude) {
        const float stepSize = 2.0f / gridSize;

        std::unordered_map<glm::ivec3, float, IVec3Comparator> grid;
        PerlinNoise perlin = PerlinNoise();

        for (int x = 0; x < gridSize; ++x) {
            for (int y = 0; y < gridSize; ++y) {
                for (int z = 0; z < gridSize; ++z) {
                    glm::vec3 position = glm::vec3(
                        (x - gridSize / 2) / float(gridSize / 2),
                        (y - gridSize / 2) / float(gridSize / 2),
                        (z - gridSize / 2) / float(gridSize / 2)
                    );

                    float baseDensity = glm::length(position);

                    float noise = 0.0f;
                    float amplitude = noiseAmplitude;
                    float frequency = noiseFrequency;

                    for (int octave = 0; octave < 6; ++octave) {
                        float ridgedNoise = 1.0f - fabs(perlin.Noise(position * frequency));
                        ridgedNoise = ridgedNoise * ridgedNoise; // Emphasize peaks
                        noise += ridgedNoise * amplitude;

                        amplitude *= 0.5f; // Reduce amplitude for finer octaves
                        frequency *= 2.0f; // Increase frequency for finer details
                    }

                    // Fade noise influence near the surface
                    float distanceToSurface = fabs(baseDensity - isoLevel);
                    float noiseFalloff = exp(-distanceToSurface * 10.0f);
                    noise *= noiseFalloff;

                    //baseDensity = -y / gridSize;
                    //baseDensity += -glm::length(position) / gridSize;
                    // Combine base sphere with noise
                    //baseDensity += perlin.Noise(glm::vec3(x, y, z));
                    //baseDensity += perlin.Noise(glm::vec3(x, y, z) * 16.12f) * 0.0575f;
                    //baseDensity += perlin.Noise(glm::vec3(x, y, z) * 8.06f) * 0.125f;
                    //baseDensity += perlin.Noise(glm::vec3(x, y, z) * 4.03f) * 0.25f;
                    //baseDensity += perlin.Noise(glm::vec3(x, y, z) * 1.96f) * 0.50f;
                    //baseDensity += perlin.Noise(glm::vec3(x, y, z) * 1.01f) * 1.00f;
                    grid[glm::vec3(x, y, z)] = baseDensity - noise;
                }
            }
        }
        return grid;
    }

    const glm::vec3 offsets[8] = {
        {-1, -1, -1},
        {1, -1, -1},
        {-1, 1, -1},
        {1, 1, -1},
        {-1, -1, 1},
        {1, -1, 1},
        {-1, 1, 1},
        {1, 1, 1}
    };

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

    void OctreeNode::Update(const glm::vec3 &playerPosition, const float planetRadius) {
        //if (IsCompletelyOutsideSphere(glm::vec3(0.0f), planetRadius)) {
        //	// Don't subdivide or update
        //	PL_INFO("Outside sphere");
        //	return;
        //}
        if (glm::length(origin) - size * 0.5f > planetRadius)
            return;

        glm::vec3 halfSize = glm::vec3(size / 2.0f);
        glm::vec3 minBounds = origin - halfSize;
        glm::vec3 maxBounds = origin + halfSize;

        bool inside =
                playerPosition.x >= minBounds.x && playerPosition.x <= maxBounds.x &&
                playerPosition.y >= minBounds.y && playerPosition.y <= maxBounds.y &&
                playerPosition.z >= minBounds.z && playerPosition.z <= maxBounds.z;

        if (size == planetRadius)
            inside = true;

        float minSize = 32.0f;
        if (inside && size > minSize) {
            if (isLeaf) SubDivide();
            for (auto &child: children) {
                if (child) child->Update(playerPosition, planetRadius);
            }
            isLeaf = false;
        } else if (!inside || size <= minSize) {
            // Too small for current distance — collapse node
            for (auto &child: children) {
                child.reset();
            }
            isLeaf = true;
        }
    }

    bool OctreeNode::IsCompletelyOutsideSphere(const glm::vec3 &center, float radius) {
        glm::vec3 halfSize = glm::vec3(size / 2.0f);
        glm::vec3 minBounds = origin - halfSize;
        glm::vec3 maxBounds = origin + halfSize;

        // Clamp center to AABB
        glm::vec3 closestPoint = glm::clamp(center, minBounds, maxBounds);
        float distance = glm::distance(center, closestPoint);

        return distance > radius;
    }
}

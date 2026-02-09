#include "pch.h"
#include "BigPhysicalBodyScript.h"
#include "PhysicalBodyScript.h"
#include "Octree.h"

namespace Plaza {
	void BigPhysicalBodyScript::OnStart(Scene* scene) {}

	void BigPhysicalBodyScript::OnUpdate(Scene* scene) {
		mPosition = scene->GetComponent<TransformComponent>(mEntityUuid)->GetWorldPosition();

		Application::Get()->mRenderer->mDebugRenderer->AddBox(
			mRootOctreeNode.mOrigin, glm::quat(), glm::vec3(mRootOctreeNode.mSize), PlColor(0, 1, 0, 0.5f), false);

		mRootOctreeNode = PlanetOctreeNode();
		mRootOctreeNode.mSize = mRadius * 2;
		mRootOctreeNode.mIsLeaf = true;
		mRootOctreeNode.mOrigin = glm::ivec3(0, 0, 0);

		mRootOctreeNode.Update(
			scene->GetComponent<TransformComponent>(CppHelper::FindEntity("PlayerEntity")->uuid)->GetWorldPosition(),
			mRadius, false);

		std::vector<PlanetOctreeNode*> renderNodes = std::vector<PlanetOctreeNode*>();
		mRootOctreeNode.GetLeafNodes(renderNodes);




		for (auto& node : renderNodes) {
			if (node == nullptr)
				continue;

			//DCOctreeNode* dcNode = static_cast<DCOctreeNode*>(node);
			//Mesh* chunkMesh = GenerateDualContouringMesh(dcNode, mIsoLevel);

			//if (chunkMesh) {
			std::string originString =
				std::to_string(node->mOrigin.x) + std::to_string(node->mOrigin.y) + std::to_string(node->mOrigin.z);

			string chunkName = "PlanetChunk " + std::to_string(node->mSize) + originString;
			Entity* foundEntity = CppHelper::FindEntity(chunkName);
			if (foundEntity != nullptr)
				scene->RemoveEntity(foundEntity->uuid);

			Entity* chunkEntity =
				scene->NewEntity(chunkName, scene->GetEntity(this->mEntityUuid));
				node->mEntityUuid = chunkEntity->uuid;
				ECS::TransformSystem::SetLocalPosition(*scene->GetComponent<TransformComponent>(chunkEntity->uuid),
													   scene,
													   glm::vec3(node->mOrigin) - glm::vec3(node->mSize * 0.5f));
				ECS::TransformSystem::SetLocalScale(*scene->GetComponent<TransformComponent>(chunkEntity->uuid), scene,
													glm::vec3(1.0f));

				//MeshRenderer* meshRenderer = scene->NewComponent<MeshRenderer>(chunkEntity->uuid);
				//meshRenderer->ChangeMesh(chunkMesh);
				//for (Material* material : mMaterialsVector) {
				//	meshRenderer->AddMaterial(material);
				//}

				// Collider* collider = scene->NewComponent<Collider>(chunkEntity->uuid);
				// collider->AddMeshShape(chunkMesh);
				// ECS::ColliderSystem::InitCollider(collider,
				// scene->GetComponent<TransformComponent>(chunkEntity->uuid), 								  nullptr);
			//}
		}


		DrawOctreeNode(mRootOctreeNode, glm::vec3(0.0f), mRadius);
	}

	void BigPhysicalBodyScript::DrawOctreeNode(const PlanetOctreeNode& node, const glm::vec3& parentOrigin,
											   float parentSize) {
		Application::Get()->mRenderer->mDebugRenderer->AddBox(node.mOrigin, glm::quat(), glm::vec3(node.mSize),
															  PlColor(0, 1, 0, 0.5f), false);
		for (const auto& child : node.mChildren) {
			if (child) {
				DrawOctreeNode(*child, node.mOrigin, node.mSize);
			}
		}
	}
} // namespace Plaza
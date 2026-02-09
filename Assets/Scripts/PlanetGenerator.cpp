#include "PlanetGenerator.h"

#include "Constants.h"
#include "Physics/PhysicalBodyScript.h"
#include "Physics/BigPhysicalBodyScript.h"
#include <thread>
#include "Engine/Core/Renderer/DebugRenderer.h"

namespace Plaza {

	void PlanetGenerator::OnStart(Scene* scene) {
		ECS::RegisterComponents();
		mMaterialsVector =
			scene->GetComponent<MeshRenderer>(CppHelper::FindEntity("MaterialsEntity")->uuid)->GetMaterials();

		GeneratePlanet(scene);
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
		scene->NewComponent<CppScriptComponent>(mPlanetEntity->uuid);
		BigPhysicalBodyScript* bigBody = static_cast<BigPhysicalBodyScript*>(
			scene->GetComponent<CppScriptComponent>(mPlanetEntity->uuid)
				->AddScriptNewInstance(scene, AssetsManager::GetScriptByName("BigPhysicalBodyScript.h")->mAssetUuid));

		mLastGeneratedPlanet = mPlanetEntity;

		double g = 9.81;
		double planetMass = g * (mPlanetRadius * mPlanetRadius) / 6.67430e-11;

		bigBody->mRadius = mPlanetRadius;
		bigBody->mWeight = planetMass;
		bigBody->mPosition = glm::vec3(0.0f);
		bigBody->mRotation = glm::quat();

		bigBody->mRootOctreeNode = PlanetOctreeNode();
		bigBody->mRootOctreeNode.mSize = mPlanetRadius * 2;
		bigBody->mRootOctreeNode.mOrigin = glm::ivec3(0, 0, 0);
		bigBody->mRootOctreeNode.mIsLeaf = true;


		//scene->GetComponent<TransformComponent>(CppHelper::FindEntity("PlayerEntity")->uuid)->GetWorldPosition()
		bigBody->mRootOctreeNode.Update(glm::vec3(0.0f, mPlanetRadius, 0.0f), mPlanetRadius, false);

		Application::Get()->mRenderer->mDebugRenderer->AddBox(bigBody->mRootOctreeNode.mOrigin, glm::quat(),
															  glm::vec3(bigBody->mRootOctreeNode.mSize),
															  PlColor(0, 1, 0, 0.5f), false);
		return mPlanetEntity;
	}

} // namespace Plaza

#include "pch.h"
#include "BigPhysicalBodyScript.h"
#include "PhysicalBodyScript.h"

namespace Plaza {
	void BigPhysicalBodyScript::OnStart(Scene* scene) {

	}

	void BigPhysicalBodyScript::OnUpdate(Scene* scene) {
		mPosition = scene->GetComponent<TransformComponent>(mEntityUuid)->GetWorldPosition();
	}
}
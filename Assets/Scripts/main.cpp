#include "pch.h"
#include "Engine/ECS/ECSManager.h"

BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved) {
    switch (ul_reason_for_call) {
    case DLL_PROCESS_ATTACH:
        PL_CORE_INFO("DLL Loaded");
        Plaza::ECS::RegisterComponents();
        break;
    case DLL_PROCESS_DETACH:
        PL_CORE_INFO("DLL Unloaded");
        break;
    }
    return TRUE;
}
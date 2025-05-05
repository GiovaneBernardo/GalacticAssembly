#pragma once


#ifdef ENGINE_EXPORT
#define PLAZA_API __declspec(dllexport)
#define GAME_API __declspec(dllimport)
#define MONO_DLL_EXPORT
#elif ENGINE_IMPORT
#define PLAZA_API  __declspec(dllimport)
#define GAME_API __declspec(dllexport)
#define MONO_DLL_IMPORT
#else
#define PLAZA_API 
#define GAME_API 
#endif

#include "Engine/Core/Engine.h"
#include "Engine/Core/PreCompiledHeaders.h"
#include "Engine/Core/Scripting/CppScriptFactory.h"
#include "Engine/Components/Core/Transform.h"
#include "Engine/Core/Scripting/CppHelper.h"
#include "Engine/Core/Time.h"
#include "Engine/Core/Input/Input.h"
#include "Engine/Core/Input/Cursor.h"

local PlazaRootDir = "C:\PlazaEngine\Engine/../"
local ProjectName = "Teste"
workspace (ProjectName)
    architecture "x64"
    startproject (ProjectName)

    configurations
    {
        "Debug",
        "Release"
    }

    flags
    {
        "MultiProcessorCompile"
    }

    outputDir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"
-- Define your project
project (ProjectName)
    kind "SharedLib" -- You can change this to WindowedApp if it's a GUI application
    language "C#"
    dotnetframework "4.7.2"

    targetdir ("Binaries")
    objdir ("Intermediates")


-- Add source files (your C# source files)
files
{
"Source/**.cs",
"Properties/**.cs"
}

-- Add references to external DLLs (change the path to your DLL)
links
{
    PlazaRootDir .. "dll/PlazaScriptCore.dll"
}

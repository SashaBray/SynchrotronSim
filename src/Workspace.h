#ifndef WORKSPACE_H
#define WORKSPACE_H
// Поиск папки workspace/ — общий для всех C++ main'ов (трекер, оптимизатор).
// Поднимается вверх от каталога .exe / argv0 / cwd, пока не найдёт подкаталог "workspace".
#include <string>
#include <vector>
#include <filesystem>
#include <system_error>
#ifdef _WIN32
#include <windows.h>
#endif

inline std::string findWorkspace(const char* argv0) {
    namespace fs = std::filesystem;
    std::error_code ec;
    std::vector<fs::path> starts;
#ifdef _WIN32
    char exePath[MAX_PATH];
    DWORD len = GetModuleFileNameA(NULL, exePath, MAX_PATH);
    if (len > 0 && len < MAX_PATH) starts.push_back(fs::path(exePath).parent_path());
#endif
    if (argv0 && *argv0) { fs::path p(argv0); if (p.has_parent_path()) starts.push_back(p.parent_path()); }
    starts.push_back(fs::current_path(ec));
    for (const auto& start : starts) {
        fs::path dir = fs::absolute(start, ec);
        for (int i = 0; i < 10; ++i) {
            fs::path cand = dir / "workspace";
            if (fs::is_directory(cand, ec)) return cand.generic_string() + "/";
            fs::path parent = dir.parent_path();
            if (parent.empty() || parent == dir) break;
            dir = parent;
        }
    }
    return "workspace/";
}

#endif

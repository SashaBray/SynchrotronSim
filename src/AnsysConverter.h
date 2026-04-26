#ifndef ANSYS_CONVERTER_H
#define ANSYS_CONVERTER_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdint>

// Структура для хранения векторов поля (B или E)
struct FieldVector {
    double x, y, z;
};

// Заголовок бинарного файла для мгновенного считывания параметров сетки
struct GridHeader {
    double minX, minY, minZ;
    double maxX, maxY, maxZ;
    double stepX, stepY, stepZ;
    int32_t nx, ny, nz;
    double nominalArg; // Ток/Напряжение, при котором снята карта
};

class AnsysConverter {
public:
    /**
     * @brief Конвертирует текстовый файл Ansys в бинарный формат .bin
     * @param rawPath Путь к исходному текстовому файлу
     * @param binPath Путь для сохранения бинарного файла
     * @param nominalArg Значение аргумента (ток/напряжение) из названия файла
     * @return true если конвертация успешна
     */
    static bool convert(const std::string& rawPath, const std::string& binPath, double nominalArg) {
        std::ifstream file(rawPath);
        if (!file.is_open()) {
            std::cerr << "[AnsysConverter] Error: Cannot open raw file " << rawPath << std::endl;
            return false;
        }

        std::string line;
        GridHeader header;
        header.nominalArg = nominalArg;

        // 1. Парсинг первой строки с параметрами сетки
        // Ожидаем формат: Grid Output Min: [0mm -100mm -100mm] Max: [1000mm 100mm 100mm] Grid Size: [5mm 5mm 5mm]
        if (!std::getline(file, line) || !parseHeader(line, header)) {
            std::cerr << "[AnsysConverter] Error: Invalid header format in " << rawPath << std::endl;
            return false;
        }

        // Вычисляем количество узлов (nx, ny, nz)
        header.nx = static_cast<int32_t>(std::round((header.maxX - header.minX) / header.stepX)) + 1;
        header.ny = static_cast<int32_t>(std::round((header.maxY - header.minY) / header.stepY)) + 1;
        header.nz = static_cast<int32_t>(std::round((header.maxZ - header.minZ) / header.stepZ)) + 1;

        // Резервируем память под плоский массив данных
        std::vector<FieldVector> data(header.nx * header.ny * header.nz, {0.0, 0.0, 0.0});

        // 2. Пропускаем строку названий колонок (X, Y, Z, Vector data...)
        std::getline(file, line);

        // 3. Построчное чтение данных
        while (std::getline(file, line)) {
            if (line.empty()) continue;

            double x, y, z, fx, fy, fz;
            // sscanf значительно быстрее потоков stringstream для больших файлов
            if (sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &fx, &fy, &fz) == 6) {
                // Переводим в метры для расчета индексов (если в файле mm)
                x /= 1000.0; y /= 1000.0; z /= 1000.0;

                // Находим индексы узла в массиве
                int i = static_cast<int>(std::round((x - header.minX) / header.stepX));
                int j = static_cast<int>(std::round((y - header.minY) / header.stepY));
                int k = static_cast<int>(std::round((z - header.minZ) / header.stepZ));

                if (i >= 0 && i < header.nx && j >= 0 && j < header.ny && k >= 0 && k < header.nz) {
                    size_t idx = static_cast<size_t>(i) * (header.ny * header.nz) +
                                 static_cast<size_t>(j) * (header.nz) +
                                 static_cast<size_t>(k);
                    data[idx] = {fx, fy, fz};
                }
            }
        }

        // 4. Запись в бинарный файл
        std::ofstream out(binPath, std::ios::binary);
        if (!out) return false;

        out.write(reinterpret_cast<const char*>(&header), sizeof(GridHeader));
        out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(FieldVector));

        std::cout << "[AnsysConverter] Successfully created: " << binPath 
                  << " (" << header.nx << "x" << header.ny << "x" << header.nz << " nodes)" << std::endl;

        return true;
    }

private:
    static bool parseHeader(const std::string& line, GridHeader& h) {
        // Извлекаем числа из строки заголовка. Предполагаем входные данные в mm.
        if (sscanf(line.c_str(), "Grid Output Min: [%lfmm %lfmm %lfmm] Max: [%lfmm %lfmm %lfmm] Grid Size: [%lfmm %lfmm %lfmm]",
            &h.minX, &h.minY, &h.minZ, &h.maxX, &h.maxY, &h.maxZ, &h.stepX, &h.stepY, &h.stepZ) == 9) {
            
            // Сразу переводим всё в метры (СИ)
            h.minX /= 1000.0; h.minY /= 1000.0; h.minZ /= 1000.0;
            h.maxX /= 1000.0; h.maxY /= 1000.0; h.maxZ /= 1000.0;
            h.stepX /= 1000.0; h.stepY /= 1000.0; h.stepZ /= 1000.0;
            return true;
        }
        return false;
    }
};

#endif // ANSYS_CONVERTER_H

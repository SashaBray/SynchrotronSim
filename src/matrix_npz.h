// matrix_npz.h — минимальный писатель numpy .npz из C++ (без внешних зависимостей).
//
// .npz = ZIP-архив (без сжатия, метод "stored") из .npy-файлов. Нужен, чтобы автономное C++ ядро писало
// траектории в формат, который читает Python-анализатор (matrix_traj_export.read_local_npz / traj_analyzer).
// Поддержаны dtype: <f8 (double), <i8 (int64), |b1 (bool), <U (unicode-скаляр) — 1-D массивы и 0-D скаляры.
#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>

namespace mc {

// CRC-32 (IEEE, как в ZIP) — таблично-независимая реализация.
inline uint32_t crc32_of(const uint8_t* data, size_t n) {
    uint32_t crc = 0xFFFFFFFFu;
    for (size_t i = 0; i < n; ++i) {
        crc ^= data[i];
        for (int k = 0; k < 8; ++k) crc = (crc >> 1) ^ (0xEDB88320u & (uint32_t)(-(int32_t)(crc & 1)));
    }
    return ~crc;
}

class NpzWriter {
public:
    // Записать один .npy-массив в архив под именем name (без расширения — .npy добавляется).
    void add(const std::string& name, const std::string& descr, const std::string& shape,
             const uint8_t* data, size_t nbytes) {
        std::string npy = buildNpy(descr, shape, data, nbytes);
        entries_.push_back({ name + ".npy", std::move(npy) });
    }
    void addVecF8(const std::string& name, const std::vector<double>& v) {
        add(name, "<f8", "(" + std::to_string(v.size()) + ",)",
            reinterpret_cast<const uint8_t*>(v.data()), v.size() * 8);
    }
    void addVecI8(const std::string& name, const std::vector<int64_t>& v) {
        add(name, "<i8", "(" + std::to_string(v.size()) + ",)",
            reinterpret_cast<const uint8_t*>(v.data()), v.size() * 8);
    }
    void addScalarF8(const std::string& name, double x) {
        add(name, "<f8", "()", reinterpret_cast<const uint8_t*>(&x), 8);
    }
    void addScalarI8(const std::string& name, int64_t x) {
        add(name, "<i8", "()", reinterpret_cast<const uint8_t*>(&x), 8);
    }
    void addScalarBool(const std::string& name, bool b) {
        uint8_t v = b ? 1 : 0; add(name, "|b1", "()", &v, 1);
    }
    void addScalarStr(const std::string& name, const std::string& s) {   // 0-D <U(len) (UCS4 LE)
        std::vector<uint8_t> buf(s.size() * 4, 0);
        for (size_t i = 0; i < s.size(); ++i) buf[i * 4] = (uint8_t)s[i];  // ASCII → кодпоинт в младшем байте
        add(name, "<U" + std::to_string(s.size()), "()", buf.data(), buf.size());
    }

    bool write(const std::string& path) const {
        std::ofstream f(path, std::ios::binary);
        if (!f) return false;
        std::vector<uint32_t> crcs, offs, sizes;
        uint32_t offset = 0;
        for (const auto& e : entries_) {                       // локальные заголовки + данные
            uint32_t crc = crc32_of(reinterpret_cast<const uint8_t*>(e.second.data()), e.second.size());
            uint32_t sz = (uint32_t)e.second.size();
            crcs.push_back(crc); offs.push_back(offset); sizes.push_back(sz);
            std::string h;
            putU32(h, 0x04034b50); putU16(h, 20); putU16(h, 0); putU16(h, 0);   // sig, ver, flags, method(stored)
            putU16(h, 0); putU16(h, 0x21);                                       // modtime, moddate (фикс.)
            putU32(h, crc); putU32(h, sz); putU32(h, sz);                        // crc, comp, uncomp
            putU16(h, (uint16_t)e.first.size()); putU16(h, 0);                   // namelen, extralen
            h += e.first;
            f.write(h.data(), h.size()); f.write(e.second.data(), e.second.size());
            offset += (uint32_t)h.size() + sz;
        }
        uint32_t cdStart = offset;
        std::string cd;
        for (size_t i = 0; i < entries_.size(); ++i) {         // центральный каталог
            putU32(cd, 0x02014b50); putU16(cd, 20); putU16(cd, 20); putU16(cd, 0); putU16(cd, 0);
            putU16(cd, 0); putU16(cd, 0x21);
            putU32(cd, crcs[i]); putU32(cd, sizes[i]); putU32(cd, sizes[i]);
            putU16(cd, (uint16_t)entries_[i].first.size());
            putU16(cd, 0); putU16(cd, 0); putU16(cd, 0); putU16(cd, 0); putU32(cd, 0);
            putU32(cd, offs[i]);
            cd += entries_[i].first;
        }
        f.write(cd.data(), cd.size());
        std::string eocd;                                       // конец центрального каталога
        putU32(eocd, 0x06054b50); putU16(eocd, 0); putU16(eocd, 0);
        putU16(eocd, (uint16_t)entries_.size()); putU16(eocd, (uint16_t)entries_.size());
        putU32(eocd, (uint32_t)cd.size()); putU32(eocd, cdStart); putU16(eocd, 0);
        f.write(eocd.data(), eocd.size());
        return (bool)f;
    }

private:
    std::vector<std::pair<std::string, std::string>> entries_;   // (имя.npy, содержимое)

    static void putU16(std::string& s, uint16_t v) { s.push_back((char)(v & 0xFF)); s.push_back((char)(v >> 8)); }
    static void putU32(std::string& s, uint32_t v) {
        for (int i = 0; i < 4; ++i) s.push_back((char)((v >> (8 * i)) & 0xFF));
    }
    static std::string buildNpy(const std::string& descr, const std::string& shape,
                                const uint8_t* data, size_t nbytes) {
        std::string hdr = "{'descr': '" + descr + "', 'fortran_order': False, 'shape': " + shape + ", }";
        size_t unpadded = 10 + hdr.size() + 1;                  // magic6+ver2+len2 + hdr + '\n'
        size_t pad = (64 - unpadded % 64) % 64;
        hdr += std::string(pad, ' '); hdr += '\n';
        std::string out;
        out += '\x93'; out += "NUMPY"; out += '\x01'; out += '\x00';   // magic + версия 1.0
        putU16(out, (uint16_t)hdr.size());
        out += hdr;
        out.append(reinterpret_cast<const char*>(data), nbytes);
        return out;
    }
};

} // namespace mc

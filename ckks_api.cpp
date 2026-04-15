#include "ckks_api.h"

#include <seal/seal.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace seal;
using namespace std;

namespace {
constexpr uint32_t kPackageVersion = 2;
constexpr uint64_t kChecksumOffset = 1469598103934665603ULL;
constexpr uint64_t kChecksumPrime = 1099511628211ULL;
constexpr int kLegacyOutputDecimals = 11;
constexpr char kPackageMagic[8] = {'C', 'K', 'K', 'S', 'P', 'K', 'G', '2'};

struct ParsedTable {
    vector<vector<double>> rows;
    vector<int> observed_column_decimals;
    size_t row_count = 0;
    size_t column_count = 0;
    size_t value_count = 0;
};

struct PackedData {
    vector<vector<double>> plaintext_blocks;
    size_t block_rows = 1;
    size_t ciphertext_count = 0;
    bool structure_pack_used = true;
};

struct EncryptedPackage {
    ckks_cipher_info info{};
    vector<int> column_decimals;
    vector<Ciphertext> ciphertexts;
};

void set_error(const string &msg, char *error_message, size_t error_message_size) {
    if (!error_message || error_message_size == 0) {
        return;
    }
    const size_t n = min(error_message_size - 1, msg.size());
    for (size_t i = 0; i < n; ++i) {
        error_message[i] = msg[i];
    }
    error_message[n] = '\0';
}

int clamp_int(int value, int low, int high) {
    return max(low, min(value, high));
}

int trim_trailing_zero_decimals(const string &token) {
    const size_t dot_pos = token.find('.');
    if (dot_pos == string::npos) {
        return 0;
    }

    const size_t exp_pos = token.find_first_of("eE");
    const size_t frac_end = (exp_pos == string::npos) ? token.size() : exp_pos;
    string fractional = token.substr(dot_pos + 1, frac_end - dot_pos - 1);
    while (!fractional.empty() && fractional.back() == '0') {
        fractional.pop_back();
    }
    return static_cast<int>(fractional.size());
}

vector<string> tokenize_line(const string &line) {
    string normalized = line;
    for (char &ch : normalized) {
        if (ch == '\t' || ch == ',') {
            ch = ' ';
        }
    }

    istringstream iss(normalized);
    vector<string> tokens;
    string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

ParsedTable read_plain_table(const string &path) {
    ifstream in(path, ios::in);
    if (!in.is_open()) {
        throw runtime_error("无法打开输入文件: " + path);
    }

    ParsedTable table;
    string line;
    while (getline(in, line)) {
        vector<string> tokens = tokenize_line(line);
        if (tokens.empty()) {
            continue;
        }

        if (table.column_count == 0) {
            table.column_count = tokens.size();
            table.observed_column_decimals.assign(table.column_count, 0);
        } else if (tokens.size() != table.column_count) {
            throw runtime_error("输入数据列数不一致，无法进行结构化打包");
        }

        vector<double> row;
        row.reserve(tokens.size());
        for (size_t i = 0; i < tokens.size(); ++i) {
            row.push_back(stod(tokens[i]));
            table.observed_column_decimals[i] = max(table.observed_column_decimals[i], trim_trailing_zero_decimals(tokens[i]));
        }
        table.rows.push_back(move(row));
    }

    if (table.rows.empty()) {
        throw runtime_error("输入文件没有可用数字: " + path);
    }

    table.row_count = table.rows.size();
    table.value_count = table.row_count * table.column_count;
    return table;
}

vector<int> choose_column_decimals(const ParsedTable &table, const ckks_encrypt_config &config) {
    vector<int> chosen(table.column_count, 0);
    const int min_decimals = max(0, config.min_quantize_decimals);
    const int default_max = (config.max_quantize_decimals > 0) ? config.max_quantize_decimals : 11;
    const int max_decimals = max(min_decimals, default_max);

    for (size_t i = 0; i < table.column_count; ++i) {
        int value = config.enable_adaptive_quantization ? table.observed_column_decimals[i] : max_decimals;
        chosen[i] = clamp_int(value, min_decimals, max_decimals);
    }
    return chosen;
}

void quantize_table_inplace(ParsedTable &table, const vector<int> &column_decimals) {
    for (size_t row = 0; row < table.row_count; ++row) {
        for (size_t col = 0; col < table.column_count; ++col) {
            const double factor = pow(10.0, static_cast<double>(column_decimals[col]));
            table.rows[row][col] = round(table.rows[row][col] * factor) / factor;
        }
    }
}

int derive_scale_bits(int max_quantize_decimals) {
    if (max_quantize_decimals <= 2) {
        return 35;
    }
    if (max_quantize_decimals <= 4) {
        return 40;
    }
    if (max_quantize_decimals <= 6) {
        return 45;
    }
    if (max_quantize_decimals <= 10) {
        return 50;
    }
    return 55;
}

EncryptionParameters build_auto_params(const ckks_encrypt_config &config) {
    const size_t expected_columns = (config.expected_columns > 0) ? config.expected_columns : 128;
    const size_t expected_rows = (config.expected_rows_per_ciphertext > 0) ? config.expected_rows_per_ciphertext : 32;
    const size_t desired_slots = max(expected_columns, expected_columns * expected_rows);

    size_t poly_modulus_degree = config.preferred_poly_modulus_degree;
    if (poly_modulus_degree == 0) {
        if (desired_slots <= 2048) {
            poly_modulus_degree = 16384;
        } else if (desired_slots <= 8192) {
            poly_modulus_degree = 32768;
        } else if (desired_slots <= 16384) {
            poly_modulus_degree = 16384;
        } else {
            poly_modulus_degree = 32768;
        }
    }

    const int max_quantize_decimals = (config.max_quantize_decimals > 0) ? config.max_quantize_decimals : 11;
    const int requested_scale_bits = (config.target_scale_bits > 0) ? config.target_scale_bits : derive_scale_bits(max_quantize_decimals);

    if (config.preferred_poly_modulus_degree == 0 && requested_scale_bits >= 55 && poly_modulus_degree < 32768) {
        poly_modulus_degree = 32768;
    }

    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    if (requested_scale_bits >= 55) {
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, 58, 58, 58, 60}));
    } else if (requested_scale_bits >= 50) {
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, 56, 56, 56, 60}));
    } else if (requested_scale_bits >= 45) {
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, 52, 52, 52, 60}));
    } else {
        const int middle_prime_bits = clamp_int(requested_scale_bits + 6, 42, 52);
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, middle_prime_bits, middle_prime_bits, 60}));
    }
    return parms;
}

size_t coeff_modulus_total_bits(const EncryptionParameters &parms) {
    size_t total = 0;
    for (const auto &modulus : parms.coeff_modulus()) {
        total += modulus.bit_count();
    }
    return total;
}

int max_supported_scale_bits(const EncryptionParameters &parms) {
    int limit = numeric_limits<int>::max();
    for (const auto &modulus : parms.coeff_modulus()) {
        limit = min(limit, static_cast<int>(modulus.bit_count()) - 2);
    }
    if (limit == numeric_limits<int>::max()) {
        limit = 40;
    }
    return max(20, limit);
}

PackedData pack_table(const ParsedTable &table, size_t slot_count, const ckks_encrypt_config &config) {
    if (table.column_count == 0) {
        throw runtime_error("输入表结构无效");
    }

    PackedData packed;
    packed.structure_pack_used = (config.enable_structure_pack != 0);

    if (!packed.structure_pack_used) {
        vector<double> linear_values;
        linear_values.reserve(table.value_count);
        for (const auto &row : table.rows) {
            linear_values.insert(linear_values.end(), row.begin(), row.end());
        }

        for (size_t offset = 0; offset < linear_values.size(); offset += slot_count) {
            const size_t chunk_size = min(slot_count, linear_values.size() - offset);
            packed.plaintext_blocks.emplace_back(linear_values.begin() + static_cast<ptrdiff_t>(offset),
                                                 linear_values.begin() + static_cast<ptrdiff_t>(offset + chunk_size));
        }
        packed.block_rows = 1;
        packed.ciphertext_count = packed.plaintext_blocks.size();
        return packed;
    }

    if (table.column_count > slot_count) {
        throw runtime_error("单条记录字段数超过 CKKS 槽位数，无法进行结构感知打包");
    }

    size_t block_rows = config.preferred_block_rows;
    const size_t max_block_rows = max<size_t>(1, slot_count / table.column_count);
    if (block_rows == 0 || block_rows > max_block_rows) {
        block_rows = max_block_rows;
    }
    block_rows = max<size_t>(1, block_rows);

    const size_t block_size = block_rows * table.column_count;
    for (size_t row_offset = 0; row_offset < table.row_count; row_offset += block_rows) {
        vector<double> block(block_size, 0.0);
        const size_t rows_in_block = min(block_rows, table.row_count - row_offset);
        for (size_t r = 0; r < rows_in_block; ++r) {
            copy(table.rows[row_offset + r].begin(),
                 table.rows[row_offset + r].end(),
                 block.begin() + static_cast<ptrdiff_t>(r * table.column_count));
        }
        packed.plaintext_blocks.push_back(move(block));
    }

    packed.block_rows = block_rows;
    packed.ciphertext_count = packed.plaintext_blocks.size();
    return packed;
}

template <typename T>
void save_binary(const string &path, const T &obj) {
    ofstream out(path, ios::binary | ios::trunc);
    if (!out.is_open()) {
        throw runtime_error("无法写入文件: " + path);
    }
    obj.save(out);
}

template <typename T>
void load_binary(const string &path, T &obj, const shared_ptr<SEALContext> &context = nullptr) {
    ifstream in(path, ios::binary);
    if (!in.is_open()) {
        throw runtime_error("无法读取文件: " + path);
    }

    if constexpr (is_same<T, Ciphertext>::value || is_same<T, SecretKey>::value || is_same<T, PublicKey>::value) {
        if (!context) {
            throw runtime_error("加载对象缺少上下文: " + path);
        }
        obj.load(*context, in);
    } else {
        obj.load(in);
    }
}

void append_u32(vector<uint8_t> &buffer, uint32_t value) {
    for (int i = 0; i < 4; ++i) {
        buffer.push_back(static_cast<uint8_t>((value >> (i * 8)) & 0xFFU));
    }
}

void append_u64(vector<uint8_t> &buffer, uint64_t value) {
    for (int i = 0; i < 8; ++i) {
        buffer.push_back(static_cast<uint8_t>((value >> (i * 8)) & 0xFFULL));
    }
}

uint32_t read_u32(const vector<uint8_t> &buffer, size_t &offset) {
    if (offset + 4 > buffer.size()) {
        throw runtime_error("密文包头损坏：读取 uint32 越界");
    }
    uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
        value |= static_cast<uint32_t>(buffer[offset + i]) << (i * 8);
    }
    offset += 4;
    return value;
}

uint64_t read_u64(const vector<uint8_t> &buffer, size_t &offset) {
    if (offset + 8 > buffer.size()) {
        throw runtime_error("密文包头损坏：读取 uint64 越界");
    }
    uint64_t value = 0;
    for (int i = 0; i < 8; ++i) {
        value |= static_cast<uint64_t>(buffer[offset + i]) << (i * 8);
    }
    offset += 8;
    return value;
}

uint64_t fnv1a64(const vector<uint8_t> &buffer, size_t size) {
    uint64_t hash = kChecksumOffset;
    for (size_t i = 0; i < size; ++i) {
        hash ^= static_cast<uint64_t>(buffer[i]);
        hash *= kChecksumPrime;
    }
    return hash;
}

vector<uint8_t> serialize_ciphertext(const Ciphertext &ciphertext) {
    stringstream stream(ios::in | ios::out | ios::binary);
    ciphertext.save(stream);
    const string payload = stream.str();
    return vector<uint8_t>(payload.begin(), payload.end());
}

Ciphertext deserialize_ciphertext(const vector<uint8_t> &payload, const shared_ptr<SEALContext> &context) {
    Ciphertext ciphertext;
    string payload_string(payload.begin(), payload.end());
    stringstream stream(payload_string, ios::in | ios::binary);
    ciphertext.load(*context, stream);
    return ciphertext;
}

void save_encrypted_package(const string &path,
                            const ckks_cipher_info &info,
                            const vector<int> &column_decimals,
                            const vector<Ciphertext> &ciphertexts) {
    vector<uint8_t> bytes;
    bytes.insert(bytes.end(), begin(kPackageMagic), end(kPackageMagic));
    append_u32(bytes, kPackageVersion);

    uint32_t flags = 0;
    if (info.adaptive_quantization_used) {
        flags |= 0x1U;
    }
    if (info.structure_pack_used) {
        flags |= 0x2U;
    }
    flags |= 0x4U;
    append_u32(bytes, flags);

    append_u64(bytes, static_cast<uint64_t>(info.row_count));
    append_u64(bytes, static_cast<uint64_t>(info.column_count));
    append_u64(bytes, static_cast<uint64_t>(info.value_count));
    append_u64(bytes, static_cast<uint64_t>(info.block_rows));
    append_u64(bytes, static_cast<uint64_t>(info.slots_per_ciphertext));
    append_u64(bytes, static_cast<uint64_t>(info.ciphertext_count));
    append_u64(bytes, static_cast<uint64_t>(info.poly_modulus_degree));
    append_u64(bytes, static_cast<uint64_t>(info.coeff_modulus_total_bits));
    append_u32(bytes, static_cast<uint32_t>(info.scale_bits));
    append_u64(bytes, static_cast<uint64_t>(column_decimals.size()));
    for (int decimals : column_decimals) {
        append_u32(bytes, static_cast<uint32_t>(max(0, decimals)));
    }

    for (const auto &ciphertext : ciphertexts) {
        vector<uint8_t> payload = serialize_ciphertext(ciphertext);
        append_u64(bytes, static_cast<uint64_t>(payload.size()));
        bytes.insert(bytes.end(), payload.begin(), payload.end());
    }

    const uint64_t checksum = fnv1a64(bytes, bytes.size());
    append_u64(bytes, checksum);

    ofstream out(path, ios::binary | ios::trunc);
    if (!out.is_open()) {
        throw runtime_error("无法写入文件: " + path);
    }
    out.write(reinterpret_cast<const char *>(bytes.data()), static_cast<streamsize>(bytes.size()));
    if (!out.good()) {
        throw runtime_error("写入密文文件失败: " + path);
    }
}

EncryptedPackage load_encrypted_package(const string &path, const shared_ptr<SEALContext> &context) {
    ifstream in(path, ios::binary);
    if (!in.is_open()) {
        throw runtime_error("无法读取文件: " + path);
    }
    vector<uint8_t> bytes((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
    if (bytes.size() < sizeof(kPackageMagic) + 4 + 4 + 8) {
        throw runtime_error("密文文件过短或损坏");
    }

    size_t offset = 0;
    for (char magic_byte : kPackageMagic) {
        if (bytes[offset++] != static_cast<uint8_t>(magic_byte)) {
            throw runtime_error("密文文件魔数不匹配");
        }
    }

    uint64_t stored_checksum = 0;
    const size_t checksum_base = bytes.size() - 8;
    for (int i = 0; i < 8; ++i) {
        stored_checksum |= static_cast<uint64_t>(bytes[checksum_base + i]) << (i * 8);
    }
    const uint64_t actual_checksum = fnv1a64(bytes, bytes.size() - 8);
    if (stored_checksum != actual_checksum) {
        throw runtime_error("密文完整性校验失败，文件可能已损坏或被篡改");
    }

    const uint32_t version = read_u32(bytes, offset);
    if (version != kPackageVersion) {
        throw runtime_error("不支持的密文包版本");
    }

    const uint32_t flags = read_u32(bytes, offset);
    EncryptedPackage package;
    package.info.row_count = static_cast<size_t>(read_u64(bytes, offset));
    package.info.column_count = static_cast<size_t>(read_u64(bytes, offset));
    package.info.value_count = static_cast<size_t>(read_u64(bytes, offset));
    package.info.block_rows = static_cast<size_t>(read_u64(bytes, offset));
    package.info.slots_per_ciphertext = static_cast<size_t>(read_u64(bytes, offset));
    package.info.ciphertext_count = static_cast<size_t>(read_u64(bytes, offset));
    package.info.poly_modulus_degree = static_cast<size_t>(read_u64(bytes, offset));
    package.info.coeff_modulus_total_bits = static_cast<size_t>(read_u64(bytes, offset));
    package.info.scale_bits = static_cast<int>(read_u32(bytes, offset));
    package.info.adaptive_quantization_used = (flags & 0x1U) ? 1 : 0;
    package.info.structure_pack_used = (flags & 0x2U) ? 1 : 0;

    const size_t decimal_count = static_cast<size_t>(read_u64(bytes, offset));
    package.column_decimals.reserve(decimal_count);
    for (size_t i = 0; i < decimal_count; ++i) {
        package.column_decimals.push_back(static_cast<int>(read_u32(bytes, offset)));
    }

    if (package.info.poly_modulus_degree != context->key_context_data()->parms().poly_modulus_degree()) {
        throw runtime_error("密文包与当前参数文件不匹配：poly_modulus_degree 不一致");
    }

    for (size_t i = 0; i < package.info.ciphertext_count; ++i) {
        const size_t payload_size = static_cast<size_t>(read_u64(bytes, offset));
        if (offset + payload_size > bytes.size() - 8) {
            throw runtime_error("密文包损坏：密文负载长度越界");
        }
        vector<uint8_t> payload(bytes.begin() + static_cast<ptrdiff_t>(offset),
                                bytes.begin() + static_cast<ptrdiff_t>(offset + payload_size));
        package.ciphertexts.push_back(deserialize_ciphertext(payload, context));
        offset += payload_size;
    }

    if (offset != bytes.size() - 8) {
        throw runtime_error("密文包存在未识别尾部数据");
    }

    return package;
}

void save_plain_table(const string &path,
                      const vector<vector<double>> &rows,
                      const vector<int> &column_decimals) {
    ofstream out(path, ios::out | ios::trunc);
    if (!out.is_open()) {
        throw runtime_error("无法写入输出文件: " + path);
    }

    for (size_t row = 0; row < rows.size(); ++row) {
        for (size_t col = 0; col < rows[row].size(); ++col) {
            const int decimals = (col < column_decimals.size()) ? column_decimals[col] : kLegacyOutputDecimals;
            const double factor = pow(10.0, static_cast<double>(decimals));
            const double rounded = (factor > 0.0) ? round(rows[row][col] * factor) / factor : rows[row][col];
            
            if (decimals == 0) {
                // 整数输出：不带小数点
                out << static_cast<long long>(round(rounded));
            } else {
                // 自适应小数位输出：只输出必要的有效位数
                ostringstream oss;
                oss << fixed << setprecision(decimals) << rounded;
                string formatted = oss.str();
                // 移除尾部多余的零，但保留至少一位小数（如果原始精度 > 0）
                size_t dot_pos = formatted.find('.');
                if (dot_pos != string::npos) {
                    size_t last_nonzero = formatted.find_last_not_of('0');
                    if (last_nonzero != string::npos && last_nonzero > dot_pos) {
                        formatted = formatted.substr(0, last_nonzero + 1);
                    } else if (last_nonzero == dot_pos) {
                        // 小数部分全为零，移除小数点
                        formatted = formatted.substr(0, dot_pos);
                    }
                }
                out << formatted;
            }
            if (col + 1 < rows[row].size()) {
                out << '\t';
            }
        }
        out << '\n';
    }
}

vector<vector<double>> unpack_blocks(const vector<vector<double>> &decoded_blocks,
                                     size_t row_count,
                                     size_t column_count,
                                     size_t block_rows,
                                     bool structure_pack_used) {
    vector<vector<double>> rows;
    if (!structure_pack_used) {
        vector<double> values;
        for (const auto &block : decoded_blocks) {
            values.insert(values.end(), block.begin(), block.end());
        }
        values.resize(row_count * column_count);
        rows.assign(row_count, vector<double>(column_count, 0.0));
        for (size_t row = 0; row < row_count; ++row) {
            copy(values.begin() + static_cast<ptrdiff_t>(row * column_count),
                 values.begin() + static_cast<ptrdiff_t>((row + 1) * column_count),
                 rows[row].begin());
        }
        return rows;
    }

    rows.reserve(row_count);
    for (size_t block_index = 0; block_index < decoded_blocks.size() && rows.size() < row_count; ++block_index) {
        for (size_t local_row = 0; local_row < block_rows && rows.size() < row_count; ++local_row) {
            const size_t offset = local_row * column_count;
            rows.emplace_back(decoded_blocks[block_index].begin() + static_cast<ptrdiff_t>(offset),
                              decoded_blocks[block_index].begin() + static_cast<ptrdiff_t>(offset + column_count));
        }
    }
    return rows;
}

ckks_encrypt_config make_default_config() {
    ckks_encrypt_config config{};
    config.expected_columns = 128;
    config.expected_rows_per_ciphertext = 32;
    config.preferred_poly_modulus_degree = 0;
    config.preferred_block_rows = 0;
    config.enable_adaptive_quantization = 1;
    config.enable_structure_pack = 1;
    config.min_quantize_decimals = 0;
    config.max_quantize_decimals = 11;
    config.target_scale_bits = 0;
    return config;
}
}  // namespace

extern "C" void ckks_init_default_config(ckks_encrypt_config *config) {
    if (!config) {
        return;
    }
    *config = make_default_config();
}

extern "C" int ckks_center_keygen_ex(const char *params_file,
                                      const char *public_key_file,
                                      const char *secret_key_file,
                                      const ckks_encrypt_config *config,
                                      char *error_message,
                                      size_t error_message_size) {
    try {
        const ckks_encrypt_config effective_config = config ? *config : make_default_config();
        EncryptionParameters parms = build_auto_params(effective_config);
        auto context = make_shared<SEALContext>(parms);
        if (!context->parameters_set()) {
            throw runtime_error("CKKS 参数无效");
        }

        KeyGenerator keygen(*context);
        PublicKey public_key;
        keygen.create_public_key(public_key);
        SecretKey secret_key = keygen.secret_key();

        save_binary(params_file, parms);
        save_binary(public_key_file, public_key);
        save_binary(secret_key_file, secret_key);
        return 0;
    } catch (const exception &e) {
        set_error(e.what(), error_message, error_message_size);
        return 1;
    }
}

extern "C" int ckks_center_keygen(const char *params_file,
                                   const char *public_key_file,
                                   const char *secret_key_file,
                                   char *error_message,
                                   size_t error_message_size) {
    return ckks_center_keygen_ex(params_file,
                                 public_key_file,
                                 secret_key_file,
                                 nullptr,
                                 error_message,
                                 error_message_size);
}

extern "C" int ckks_user_encrypt_ex(const char *input_file,
                                     const char *params_file,
                                     const char *public_key_file,
                                     const char *encrypted_file,
                                     const ckks_encrypt_config *config,
                                     ckks_cipher_info *cipher_info,
                                     char *error_message,
                                     size_t error_message_size) {
    try {
        const ckks_encrypt_config effective_config = config ? *config : make_default_config();
        ParsedTable table = read_plain_table(input_file);
        vector<int> column_decimals = choose_column_decimals(table, effective_config);
        quantize_table_inplace(table, column_decimals);

        EncryptionParameters parms;
        load_binary(params_file, parms);
        auto context = make_shared<SEALContext>(parms);
        if (!context->parameters_set()) {
            throw runtime_error("加载后的 CKKS 参数无效");
        }

        PublicKey public_key;
        load_binary(public_key_file, public_key, context);

        CKKSEncoder encoder(*context);
        Encryptor encryptor(*context, public_key);

        const size_t slot_count = encoder.slot_count();
        PackedData packed = pack_table(table, slot_count, effective_config);
        const int max_decimals = column_decimals.empty() ? 0 : *max_element(column_decimals.begin(), column_decimals.end());
        const int desired_scale_bits = (effective_config.target_scale_bits > 0)
                                           ? effective_config.target_scale_bits
                                           : derive_scale_bits(max_decimals);
        const int scale_bits = min(desired_scale_bits, max_supported_scale_bits(parms));
        const double scale = pow(2.0, static_cast<double>(scale_bits));

        vector<Ciphertext> ciphertexts;
        ciphertexts.reserve(packed.plaintext_blocks.size());
        for (const auto &block : packed.plaintext_blocks) {
            Plaintext plain;
            encoder.encode(block, scale, plain);
            Ciphertext encrypted;
            encryptor.encrypt(plain, encrypted);
            ciphertexts.push_back(move(encrypted));
        }

        ckks_cipher_info info{};
        info.row_count = table.row_count;
        info.column_count = table.column_count;
        info.value_count = table.value_count;
        info.block_rows = packed.block_rows;
        info.slots_per_ciphertext = slot_count;
        info.ciphertext_count = ciphertexts.size();
        info.poly_modulus_degree = parms.poly_modulus_degree();
        info.coeff_modulus_total_bits = coeff_modulus_total_bits(parms);
        info.scale_bits = scale_bits;
        info.adaptive_quantization_used = effective_config.enable_adaptive_quantization ? 1 : 0;
        info.structure_pack_used = packed.structure_pack_used ? 1 : 0;

        save_encrypted_package(encrypted_file, info, column_decimals, ciphertexts);
        if (cipher_info) {
            *cipher_info = info;
        }
        return 0;
    } catch (const exception &e) {
        set_error(e.what(), error_message, error_message_size);
        return 1;
    }
}

extern "C" int ckks_user_encrypt(const char *input_file,
                                  const char *params_file,
                                  const char *public_key_file,
                                  const char *encrypted_file,
                                  char *error_message,
                                  size_t error_message_size) {
    return ckks_user_encrypt_ex(input_file,
                                params_file,
                                public_key_file,
                                encrypted_file,
                                nullptr,
                                nullptr,
                                error_message,
                                error_message_size);
}

extern "C" int ckks_center_decrypt_ex(const char *params_file,
                                       const char *secret_key_file,
                                       const char *encrypted_file,
                                       const char *output_file,
                                       ckks_cipher_info *cipher_info,
                                       char *error_message,
                                       size_t error_message_size) {
    try {
        EncryptionParameters parms;
        load_binary(params_file, parms);

        auto context = make_shared<SEALContext>(parms);
        if (!context->parameters_set()) {
            throw runtime_error("加载后的 CKKS 参数无效");
        }

        SecretKey secret_key;
        load_binary(secret_key_file, secret_key, context);
        EncryptedPackage package = load_encrypted_package(encrypted_file, context);

        Decryptor decryptor(*context, secret_key);
        CKKSEncoder encoder(*context);
        vector<vector<double>> decoded_blocks;
        decoded_blocks.reserve(package.ciphertexts.size());

        for (const auto &ciphertext : package.ciphertexts) {
            Plaintext plain;
            decryptor.decrypt(ciphertext, plain);

            vector<double> values;
            encoder.decode(plain, values);
            if (package.info.structure_pack_used) {
                values.resize(package.info.block_rows * package.info.column_count);
            }
            decoded_blocks.push_back(move(values));
        }

        vector<vector<double>> rows = unpack_blocks(decoded_blocks,
                                                    package.info.row_count,
                                                    package.info.column_count,
                                                    package.info.block_rows,
                                                    package.info.structure_pack_used != 0);
        save_plain_table(output_file, rows, package.column_decimals);
        if (cipher_info) {
            *cipher_info = package.info;
        }
        return 0;
    } catch (const exception &e) {
        set_error(e.what(), error_message, error_message_size);
        return 1;
    }
}

extern "C" int ckks_center_decrypt(const char *params_file,
                                    const char *secret_key_file,
                                    const char *encrypted_file,
                                    const char *output_file,
                                    char *error_message,
                                    size_t error_message_size) {
    return ckks_center_decrypt_ex(params_file,
                                  secret_key_file,
                                  encrypted_file,
                                  output_file,
                                  nullptr,
                                  error_message,
                                  error_message_size);
}

/* ==================== 基准测试实现 ==================== */

namespace {

#ifdef _WIN32
#include <windows.h>
double get_time_ms() {
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return static_cast<double>(counter.QuadPart) * 1000.0 / static_cast<double>(freq.QuadPart);
}
#else
#include <sys/time.h>
double get_time_ms() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return static_cast<double>(tv.tv_sec) * 1000.0 + static_cast<double>(tv.tv_usec) / 1000.0;
}
#endif

size_t get_file_size(const string &path) {
    ifstream in(path, ios::binary | ios::ate);
    if (!in.is_open()) {
        return 0;
    }
    return static_cast<size_t>(in.tellg());
}

void compute_error_metrics(const ParsedTable &original,
                           const vector<vector<double>> &decrypted,
                           const vector<int> &column_decimals,
                           double &mse,
                           double &mae,
                           double &max_error) {
    mse = 0.0;
    mae = 0.0;
    max_error = 0.0;
    size_t count = 0;

    for (size_t row = 0; row < original.row_count && row < decrypted.size(); ++row) {
        for (size_t col = 0; col < original.column_count && col < decrypted[row].size(); ++col) {
            // 按照实际输出精度四舍五入后再比较（公平评估）
            const int decimals = (col < column_decimals.size()) ? column_decimals[col] : 11;
            const double factor = pow(10.0, static_cast<double>(decimals));
            const double rounded_decrypted = round(decrypted[row][col] * factor) / factor;
            const double rounded_original = round(original.rows[row][col] * factor) / factor;
            
            const double error = fabs(rounded_original - rounded_decrypted);
            mse += error * error;
            mae += error;
            max_error = max(max_error, error);
            ++count;
        }
    }

    if (count > 0) {
        mse /= static_cast<double>(count);
        mae /= static_cast<double>(count);
    }
}

int run_benchmark_internal(const string &input_file,
                           const ckks_encrypt_config &config,
                           ckks_benchmark_result &result) {
    memset(&result, 0, sizeof(result));

    // 读取原始数据
    ParsedTable original_table = read_plain_table(input_file);

    // 临时文件路径
    const string params_file = "benchmark_params.bin";
    const string public_key_file = "benchmark_public.bin";
    const string secret_key_file = "benchmark_secret.bin";
    const string encrypted_file = "benchmark_encrypted.bin";
    const string output_file = "benchmark_output.txt";

    // 密钥生成计时
    double t0 = get_time_ms();
    EncryptionParameters parms = build_auto_params(config);
    auto context = make_shared<SEALContext>(parms);
    if (!context->parameters_set()) {
        throw runtime_error("CKKS 参数无效");
    }
    KeyGenerator keygen(*context);
    PublicKey public_key;
    keygen.create_public_key(public_key);
    SecretKey secret_key = keygen.secret_key();
    save_binary(params_file, parms);
    save_binary(public_key_file, public_key);
    save_binary(secret_key_file, secret_key);
    double t1 = get_time_ms();
    result.keygen_time_ms = t1 - t0;

    // 加密计时
    t0 = get_time_ms();
    ParsedTable table = read_plain_table(input_file);
    vector<int> column_decimals = choose_column_decimals(table, config);
    quantize_table_inplace(table, column_decimals);

    CKKSEncoder encoder(*context);
    Encryptor encryptor(*context, public_key);
    const size_t slot_count = encoder.slot_count();
    PackedData packed = pack_table(table, slot_count, config);
    
    const int max_decimals = column_decimals.empty() ? 0 : *max_element(column_decimals.begin(), column_decimals.end());
    const int desired_scale_bits = (config.target_scale_bits > 0) ? config.target_scale_bits : derive_scale_bits(max_decimals);
    const int scale_bits = min(desired_scale_bits, max_supported_scale_bits(parms));
    const double scale = pow(2.0, static_cast<double>(scale_bits));

    vector<Ciphertext> ciphertexts;
    ciphertexts.reserve(packed.plaintext_blocks.size());
    for (const auto &block : packed.plaintext_blocks) {
        Plaintext plain;
        encoder.encode(block, scale, plain);
        Ciphertext encrypted;
        encryptor.encrypt(plain, encrypted);
        ciphertexts.push_back(move(encrypted));
    }

    ckks_cipher_info info{};
    info.row_count = table.row_count;
    info.column_count = table.column_count;
    info.value_count = table.value_count;
    info.block_rows = packed.block_rows;
    info.slots_per_ciphertext = slot_count;
    info.ciphertext_count = ciphertexts.size();
    info.poly_modulus_degree = parms.poly_modulus_degree();
    info.coeff_modulus_total_bits = coeff_modulus_total_bits(parms);
    info.scale_bits = scale_bits;
    info.adaptive_quantization_used = config.enable_adaptive_quantization ? 1 : 0;
    info.structure_pack_used = packed.structure_pack_used ? 1 : 0;

    save_encrypted_package(encrypted_file, info, column_decimals, ciphertexts);
    t1 = get_time_ms();
    result.encrypt_time_ms = t1 - t0;

    // 解密计时
    t0 = get_time_ms();
    EncryptedPackage package = load_encrypted_package(encrypted_file, context);
    Decryptor decryptor(*context, secret_key);
    vector<vector<double>> decoded_blocks;
    decoded_blocks.reserve(package.ciphertexts.size());

    for (const auto &ciphertext : package.ciphertexts) {
        Plaintext plain;
        decryptor.decrypt(ciphertext, plain);
        vector<double> values;
        encoder.decode(plain, values);
        if (package.info.structure_pack_used) {
            values.resize(package.info.block_rows * package.info.column_count);
        }
        decoded_blocks.push_back(move(values));
    }

    vector<vector<double>> decrypted_rows = unpack_blocks(decoded_blocks,
                                                          package.info.row_count,
                                                          package.info.column_count,
                                                          package.info.block_rows,
                                                          package.info.structure_pack_used != 0);
    save_plain_table(output_file, decrypted_rows, package.column_decimals);
    t1 = get_time_ms();
    result.decrypt_time_ms = t1 - t0;

    // 计算误差指标（按实际输出精度四舍五入后比较，公平评估）
    compute_error_metrics(original_table, decrypted_rows, column_decimals, result.mse, result.mae, result.max_error);

    // 文件大小
    result.encrypted_file_size = get_file_size(encrypted_file);
    result.params_file_size = get_file_size(params_file);
    result.public_key_file_size = get_file_size(public_key_file);
    result.total_time_ms = result.keygen_time_ms + result.encrypt_time_ms + result.decrypt_time_ms;
    result.cipher_info = info;

    // 清理临时文件
    remove(params_file.c_str());
    remove(public_key_file.c_str());
    remove(secret_key_file.c_str());
    remove(encrypted_file.c_str());
    remove(output_file.c_str());

    return 0;
}

}  // namespace

extern "C" int ckks_run_benchmark(const char *input_file,
                                   const ckks_encrypt_config *config,
                                   ckks_benchmark_result *result,
                                   char *error_message,
                                   size_t error_message_size) {
    try {
        if (!result) {
            throw runtime_error("result 指针为空");
        }
        const ckks_encrypt_config effective_config = config ? *config : make_default_config();
        return run_benchmark_internal(input_file, effective_config, *result);
    } catch (const exception &e) {
        set_error(e.what(), error_message, error_message_size);
        return 1;
    }
}

extern "C" int ckks_run_comparison(const char *input_file,
                                    ckks_comparison_result *result,
                                    char *error_message,
                                    size_t error_message_size) {
    try {
        if (!result) {
            throw runtime_error("result 指针为空");
        }
        memset(result, 0, sizeof(*result));

        // 传统模式配置：固定精度、线性打包
        ckks_encrypt_config traditional_config = make_default_config();
        traditional_config.enable_adaptive_quantization = 0;  // 关闭自适应量化
        traditional_config.enable_structure_pack = 0;         // 关闭结构化打包
        traditional_config.max_quantize_decimals = 11;        // 固定 11 位小数

        // 创新模式配置：自适应量化、结构化打包
        ckks_encrypt_config innovative_config = make_default_config();
        innovative_config.enable_adaptive_quantization = 1;   // 启用自适应量化
        innovative_config.enable_structure_pack = 1;          // 启用结构化打包

        // 运行传统模式基准测试
        int rc = run_benchmark_internal(input_file, traditional_config, result->traditional);
        if (rc != 0) {
            throw runtime_error("传统模式基准测试失败");
        }

        // 运行创新模式基准测试
        rc = run_benchmark_internal(input_file, innovative_config, result->innovative);
        if (rc != 0) {
            throw runtime_error("创新模式基准测试失败");
        }

        // 计算对比指标
        if (result->traditional.mae > 0) {
            result->precision_improvement_pct = 
                (result->traditional.mae - result->innovative.mae) / result->traditional.mae * 100.0;
        }
        if (result->traditional.encrypted_file_size > 0) {
            result->memory_reduction_pct = 
                (static_cast<double>(result->traditional.encrypted_file_size) - 
                 static_cast<double>(result->innovative.encrypted_file_size)) /
                static_cast<double>(result->traditional.encrypted_file_size) * 100.0;
        }
        if (result->traditional.total_time_ms > 0) {
            result->time_reduction_pct = 
                (result->traditional.total_time_ms - result->innovative.total_time_ms) /
                result->traditional.total_time_ms * 100.0;
        }

        return 0;
    } catch (const exception &e) {
        set_error(e.what(), error_message, error_message_size);
        return 1;
    }
}

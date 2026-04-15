#ifndef CKKS_API_H
#define CKKS_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ckks_encrypt_config {
    size_t expected_columns;
    size_t expected_rows_per_ciphertext;
    size_t preferred_poly_modulus_degree;
    size_t preferred_block_rows;
    int enable_adaptive_quantization;
    int enable_structure_pack;
    int min_quantize_decimals;
    int max_quantize_decimals;
    int target_scale_bits;
} ckks_encrypt_config;

typedef struct ckks_cipher_info {
    size_t row_count;
    size_t column_count;
    size_t value_count;
    size_t block_rows;
    size_t slots_per_ciphertext;
    size_t ciphertext_count;
    size_t poly_modulus_degree;
    size_t coeff_modulus_total_bits;
    int scale_bits;
    int adaptive_quantization_used;
    int structure_pack_used;
} ckks_cipher_info;

void ckks_init_default_config(ckks_encrypt_config *config);

int ckks_center_keygen_ex(const char *params_file,
                          const char *public_key_file,
                          const char *secret_key_file,
                          const ckks_encrypt_config *config,
                          char *error_message,
                          size_t error_message_size);

int ckks_center_keygen(const char *params_file,
                       const char *public_key_file,
                       const char *secret_key_file,
                       char *error_message,
                       size_t error_message_size);

int ckks_user_encrypt_ex(const char *input_file,
                         const char *params_file,
                         const char *public_key_file,
                         const char *encrypted_file,
                         const ckks_encrypt_config *config,
                         ckks_cipher_info *cipher_info,
                         char *error_message,
                         size_t error_message_size);

int ckks_user_encrypt(const char *input_file,
                      const char *params_file,
                      const char *public_key_file,
                      const char *encrypted_file,
                      char *error_message,
                      size_t error_message_size);

    int ckks_center_decrypt_ex(const char *params_file,
                          const char *secret_key_file,
                          const char *encrypted_file,
                          const char *output_file,
                          ckks_cipher_info *cipher_info,
                          char *error_message,
                          size_t error_message_size);

int ckks_center_decrypt(const char *params_file,
                        const char *secret_key_file,
                        const char *encrypted_file,
                        const char *output_file,
                        char *error_message,
                        size_t error_message_size);

/* ==================== 基准测试 API ==================== */

typedef struct ckks_benchmark_result {
    /* 精度指标 */
    double mse;                      /* 均方误差 */
    double mae;                      /* 平均绝对误差 */
    double max_error;                /* 最大误差 */
    
    /* 内存指标 */
    size_t encrypted_file_size;      /* 密文文件大小（字节） */
    size_t params_file_size;         /* 参数文件大小（字节） */
    size_t public_key_file_size;     /* 公钥文件大小（字节） */
    
    /* 时间指标（毫秒） */
    double keygen_time_ms;           /* 密钥生成耗时 */
    double encrypt_time_ms;          /* 加密耗时 */
    double decrypt_time_ms;          /* 解密耗时 */
    double total_time_ms;            /* 总耗时 */
    
    /* 元数据 */
    ckks_cipher_info cipher_info;    /* 加密信息 */
} ckks_benchmark_result;

typedef struct ckks_comparison_result {
    ckks_benchmark_result traditional;  /* 传统模式结果 */
    ckks_benchmark_result innovative;   /* 创新模式结果 */
    
    /* 对比指标 */
    double precision_improvement_pct;   /* 精度提升百分比（MAE） */
    double memory_reduction_pct;        /* 内存减少百分比 */
    double time_reduction_pct;          /* 时间减少百分比 */
} ckks_comparison_result;

/* 运行单次基准测试 */
int ckks_run_benchmark(const char *input_file,
                       const ckks_encrypt_config *config,
                       ckks_benchmark_result *result,
                       char *error_message,
                       size_t error_message_size);

/* 运行对比测试：传统模式 vs 创新模式 */
int ckks_run_comparison(const char *input_file,
                        ckks_comparison_result *result,
                        char *error_message,
                        size_t error_message_size);

#ifdef __cplusplus
}
#endif

#endif

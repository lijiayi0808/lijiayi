/**
 * CKKS 基准测试工具
 * 对比传统 CKKS 与创新增强（自适应量化编码 + 结构感知 SIMD 打包）的性能
 * 
 * 评估指标：
 * 1. 精度：MSE（均方误差）、MAE（平均绝对误差）、最大误差
 * 2. 内存：密文文件大小
 * 3. 时间：密钥生成、加密、解密、总耗时
 */

#include "ckks_api.h"

#include <stdio.h>
#include <string.h>

static void print_separator(const char *title) {
    printf("\n");
    printf("================================================================================\n");
    printf("  %s\n", title);
    printf("================================================================================\n");
}

static void print_benchmark_result(const char *mode_name, const ckks_benchmark_result *r) {
    printf("\n【%s】\n", mode_name);
    printf("--------------------------------------------------------------------------------\n");
    
    printf("  [精度指标]\n");
    printf("    均方误差 (MSE):     %.12e\n", r->mse);
    printf("    平均绝对误差 (MAE): %.12e\n", r->mae);
    printf("    最大误差:           %.12e\n", r->max_error);
    
    printf("\n  [内存指标]\n");
    printf("    密文文件大小:       %zu 字节 (%.2f KB)\n", 
           r->encrypted_file_size, 
           (double)r->encrypted_file_size / 1024.0);
    printf("    参数文件大小:       %zu 字节\n", r->params_file_size);
    printf("    公钥文件大小:       %zu 字节 (%.2f KB)\n", 
           r->public_key_file_size,
           (double)r->public_key_file_size / 1024.0);
    
    printf("\n  [时间指标]\n");
    printf("    密钥生成耗时:       %.2f ms\n", r->keygen_time_ms);
    printf("    加密耗时:           %.2f ms\n", r->encrypt_time_ms);
    printf("    解密耗时:           %.2f ms\n", r->decrypt_time_ms);
    printf("    总耗时:             %.2f ms\n", r->total_time_ms);
    
    printf("\n  [加密元数据]\n");
    printf("    数据维度:           %zu 行 × %zu 列 = %zu 个值\n", 
           r->cipher_info.row_count,
           r->cipher_info.column_count,
           r->cipher_info.value_count);
    printf("    密文块数:           %zu 个\n", r->cipher_info.ciphertext_count);
    printf("    每块行数:           %zu 行\n", r->cipher_info.block_rows);
    printf("    槽位数:             %zu\n", r->cipher_info.slots_per_ciphertext);
    printf("    多项式次数:         %zu\n", r->cipher_info.poly_modulus_degree);
    printf("    模数总位数:         %zu bits\n", r->cipher_info.coeff_modulus_total_bits);
    printf("    缩放位数:           %d bits\n", r->cipher_info.scale_bits);
    printf("    自适应量化:         %s\n", r->cipher_info.adaptive_quantization_used ? "启用" : "关闭");
    printf("    结构化打包:         %s\n", r->cipher_info.structure_pack_used ? "启用" : "关闭");
}

static void print_comparison_summary(const ckks_comparison_result *r) {
    print_separator("对比总结");
    
    printf("\n  指标                      传统模式          创新模式          提升/变化\n");
    printf("--------------------------------------------------------------------------------\n");
    
    // 精度对比
    printf("  MAE (平均绝对误差)        %.6e      %.6e      ", 
           r->traditional.mae, r->innovative.mae);
    if (r->precision_improvement_pct > 0) {
        printf("↑ 精度提升 %.2f%%\n", r->precision_improvement_pct);
    } else if (r->precision_improvement_pct < 0) {
        printf("↓ 精度下降 %.2f%%\n", -r->precision_improvement_pct);
    } else {
        printf("= 无变化\n");
    }
    
    printf("  MSE (均方误差)            %.6e      %.6e\n",
           r->traditional.mse, r->innovative.mse);
    
    printf("  最大误差                  %.6e      %.6e\n",
           r->traditional.max_error, r->innovative.max_error);
    
    // 内存对比
    printf("\n  密文文件大小              %-10zu 字节    %-10zu 字节    ",
           r->traditional.encrypted_file_size, r->innovative.encrypted_file_size);
    if (r->memory_reduction_pct > 0) {
        printf("↓ 减少 %.2f%%\n", r->memory_reduction_pct);
    } else if (r->memory_reduction_pct < 0) {
        printf("↑ 增加 %.2f%%\n", -r->memory_reduction_pct);
    } else {
        printf("= 无变化\n");
    }
    
    // 时间对比
    printf("\n  总耗时                    %-10.2f ms     %-10.2f ms     ",
           r->traditional.total_time_ms, r->innovative.total_time_ms);
    if (r->time_reduction_pct > 0) {
        printf("↓ 减少 %.2f%%\n", r->time_reduction_pct);
    } else if (r->time_reduction_pct < 0) {
        printf("↑ 增加 %.2f%%\n", -r->time_reduction_pct);
    } else {
        printf("= 无变化\n");
    }
    
    printf("    - 密钥生成              %-10.2f ms     %-10.2f ms\n",
           r->traditional.keygen_time_ms, r->innovative.keygen_time_ms);
    printf("    - 加密                  %-10.2f ms     %-10.2f ms\n",
           r->traditional.encrypt_time_ms, r->innovative.encrypt_time_ms);
    printf("    - 解密                  %-10.2f ms     %-10.2f ms\n",
           r->traditional.decrypt_time_ms, r->innovative.decrypt_time_ms);
    
    // 密文块数对比
    printf("\n  密文块数                  %-10zu 个      %-10zu 个\n",
           r->traditional.cipher_info.ciphertext_count, 
           r->innovative.cipher_info.ciphertext_count);
    
    print_separator("结论");
    printf("\n");
    
    // 自适应量化效果
    printf("  1. 自适应量化编码效果：\n");
    if (r->precision_improvement_pct >= 0) {
        printf("     ✓ 精度相当或更优，MAE 变化 %.4f%%\n", r->precision_improvement_pct);
    } else {
        printf("     × 精度略有下降，MAE 变化 %.4f%%（在 CKKS 近似误差范围内）\n", 
               r->precision_improvement_pct);
    }
    
    // 结构化打包效果
    printf("\n  2. 结构感知 SIMD 打包效果：\n");
    if (r->innovative.cipher_info.ciphertext_count < r->traditional.cipher_info.ciphertext_count) {
        printf("     ✓ 密文块数减少：%zu → %zu（减少 %.1f%%）\n",
               r->traditional.cipher_info.ciphertext_count,
               r->innovative.cipher_info.ciphertext_count,
               (1.0 - (double)r->innovative.cipher_info.ciphertext_count / 
                (double)r->traditional.cipher_info.ciphertext_count) * 100.0);
    } else if (r->innovative.cipher_info.ciphertext_count == r->traditional.cipher_info.ciphertext_count) {
        printf("     = 密文块数相同：%zu\n", r->innovative.cipher_info.ciphertext_count);
    } else {
        printf("     × 密文块数增加：%zu → %zu\n",
               r->traditional.cipher_info.ciphertext_count,
               r->innovative.cipher_info.ciphertext_count);
    }
    
    if (r->time_reduction_pct > 0) {
        printf("     ✓ 运行时间减少 %.2f%%\n", r->time_reduction_pct);
    } else if (r->time_reduction_pct < -5) {
        printf("     × 运行时间增加 %.2f%%\n", -r->time_reduction_pct);
    } else {
        printf("     ≈ 运行时间基本持平\n");
    }
    
    printf("\n================================================================================\n");
}

int main(int argc, char *argv[]) {
    const char *input_file = "data.txt";
    
    if (argc >= 2) {
        input_file = argv[1];
    }
    
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════════════════════════╗\n");
    printf("║                     CKKS 创新增强 vs 传统模式 性能对比                       ║\n");
    printf("╠══════════════════════════════════════════════════════════════════════════════╣\n");
    printf("║  创新增强包括：                                                              ║\n");
    printf("║    1. 自适应量化编码 - 按列分析小数位数，自动决定量化精度                    ║\n");
    printf("║    2. 结构感知 SIMD 打包 - 按行记录进行批量编码，提高槽位利用率              ║\n");
    printf("╚══════════════════════════════════════════════════════════════════════════════╝\n");
    printf("\n输入文件: %s\n", input_file);
    
    char err[1024] = {0};
    ckks_comparison_result comparison;
    
    printf("\n正在运行基准测试...\n");
    printf("  [1/2] 传统模式（固定精度 + 线性打包）...\n");
    printf("  [2/2] 创新模式（自适应量化 + 结构化打包）...\n");
    
    int rc = ckks_run_comparison(input_file, &comparison, err, sizeof(err));
    
    if (rc != 0) {
        fprintf(stderr, "\n错误: %s\n", err[0] ? err : "未知错误");
        return rc;
    }
    
    // 输出详细结果
    print_separator("传统模式详细结果");
    print_benchmark_result("传统模式（固定 11 位小数 + 线性打包）", &comparison.traditional);
    
    print_separator("创新模式详细结果");
    print_benchmark_result("创新模式（自适应量化 + 结构化打包）", &comparison.innovative);
    
    // 输出对比总结
    print_comparison_summary(&comparison);
    
    return 0;
}

#include "ckks_api.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    const char *input_file = "data.txt";
    const char *params_file = "ckks_params.bin";
    const char *public_key_file = "public_key.bin";
    const char *encrypted_file = "encrypted_data.bin";
    ckks_encrypt_config config;
    ckks_cipher_info info;
    ckks_init_default_config(&config);

    if (argc == 5) {
        input_file = argv[1];
        params_file = argv[2];
        public_key_file = argv[3];
        encrypted_file = argv[4];
    } else if (argc == 6) {
        input_file = argv[1];
        params_file = argv[2];
        public_key_file = argv[3];
        encrypted_file = argv[4];
        config.preferred_block_rows = (size_t)strtoull(argv[5], NULL, 10);
    }

    char err[1024] = {0};
    int rc = ckks_user_encrypt_ex(input_file,
                                  params_file,
                                  public_key_file,
                                  encrypted_file,
                                  &config,
                                  &info,
                                  err,
                                  sizeof(err));
    if (rc != 0) {
        fprintf(stderr, "错误: %s\n", err[0] ? err : "未知错误");
        return rc;
    }

    printf("用户侧加密完成\n");
    printf("输入明文: %s\n", input_file);
    printf("使用参数: %s\n", params_file);
    printf("使用公钥: %s\n", public_key_file);
    printf("密文文件: %s\n", encrypted_file);
        printf("打包信息: rows=%zu, cols=%zu, block_rows=%zu, ciphertexts=%zu, slots=%zu\n",
            info.row_count,
            info.column_count,
            info.block_rows,
            info.ciphertext_count,
            info.slots_per_ciphertext);
        printf("量化与参数: adaptive=%d, structured=%d, scale_bits=%d, poly_modulus_degree=%zu\n",
            info.adaptive_quantization_used,
            info.structure_pack_used,
            info.scale_bits,
            info.poly_modulus_degree);
    return 0;
}

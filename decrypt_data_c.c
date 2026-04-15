#include "ckks_api.h"

#include <stdio.h>

int main(int argc, char *argv[]) {
    const char *params_file = "ckks_params.bin";
    const char *secret_key_file = "secret_key.bin";
    const char *encrypted_file = "encrypted_data.bin";
    const char *output_file = "decrypted_data.txt";
    ckks_cipher_info info;

    if (argc == 5) {
        params_file = argv[1];
        secret_key_file = argv[2];
        encrypted_file = argv[3];
        output_file = argv[4];
    }

    char err[1024] = {0};
    int rc = ckks_center_decrypt_ex(params_file,
                                    secret_key_file,
                                    encrypted_file,
                                    output_file,
                                    &info,
                                    err,
                                    sizeof(err));
    if (rc != 0) {
        fprintf(stderr, "错误: %s\n", err[0] ? err : "未知错误");
        return rc;
    }

    printf("调度中心解密完成\n");
    printf("密文输入: %s\n", encrypted_file);
    printf("明文输出: %s\n", output_file);
    printf("密文元数据: rows=%zu, cols=%zu, ciphertexts=%zu, block_rows=%zu\n",
           info.row_count,
           info.column_count,
           info.ciphertext_count,
           info.block_rows);
    printf("参数摘要: scale_bits=%d, poly_modulus_degree=%zu, coeff_modulus_total_bits=%zu\n",
           info.scale_bits,
           info.poly_modulus_degree,
           info.coeff_modulus_total_bits);
    return 0;
}

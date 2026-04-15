#include "ckks_api.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    const char *params_file = "ckks_params.bin";
    const char *public_key_file = "public_key.bin";
    const char *secret_key_file = "secret_key.bin";
    ckks_encrypt_config config;
    ckks_init_default_config(&config);

    if (argc == 4) {
        params_file = argv[1];
        public_key_file = argv[2];
        secret_key_file = argv[3];
    } else if (argc == 7) {
        params_file = argv[1];
        public_key_file = argv[2];
        secret_key_file = argv[3];
        config.expected_columns = (size_t)strtoull(argv[4], NULL, 10);
        config.expected_rows_per_ciphertext = (size_t)strtoull(argv[5], NULL, 10);
        config.max_quantize_decimals = atoi(argv[6]);
    }

    char err[1024] = {0};
    int rc = ckks_center_keygen_ex(params_file, public_key_file, secret_key_file, &config, err, sizeof(err));
    if (rc != 0) {
        fprintf(stderr, "错误: %s\n", err[0] ? err : "未知错误");
        return rc;
    }

    printf("调度中心密钥生成完成\n");
    printf("参数文件: %s\n", params_file);
    printf("公钥文件: %s（下发给用户侧）\n", public_key_file);
    printf("私钥文件: %s（仅中心保留）\n", secret_key_file);
    printf("自动配置: expected_columns=%zu, expected_rows_per_ciphertext=%zu, max_quantize_decimals=%d\n",
           config.expected_columns,
           config.expected_rows_per_ciphertext,
           config.max_quantize_decimals);
    return 0;
}

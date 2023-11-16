#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>

#define MIN(a, b) ((a) < (b)) ? (a) : (b)

enum
{
    NO_BIT_ERROR,
    ONE_BIT_ERROR,
    TWO_BIT_ERROR,
    UNKNOWN_ERROR
};

typedef struct HAMMING_INFO_T_
{
    uint8_t **parity;
    uint32_t *hash_table;
    uint32_t n;
    uint32_t m;
    uint32_t code_len;
    uint32_t ded_on;
    uint32_t err_index;
} HAMMING_INFO_T;

typedef struct ERROR_INFO_T_
{
    uint32_t error_type;
    uint32_t error_index;
    uint32_t error_syndrome;
} ERROR_INFO_T;

uint32_t gen_random_number(void)
{
    return (rand() << 16) | (rand());
}

uint32_t is_power_of_2(uint32_t n)
{
    return n > 0 && !(n & n - 1);
}

void show_matrix(HAMMING_INFO_T *info)
{
    for (uint32_t i = 0; i < info->n; i++)
    {
        for (uint32_t j = 0; j < info->m; j++)
        {
            printf("%d, ", info->parity[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void show_table(HAMMING_INFO_T *info)
{
    uint32_t l = 1 << info->n;
    for (uint32_t i = 0; i < l; i++)
    {
        printf("%d: %d\n", i, info->hash_table[i]);
    }
}

uint32_t cal_parity_count(uint32_t len)
{
    uint32_t p_cnt = 0;
    while ((1 << p_cnt) < (len + p_cnt + 1))
    {
        p_cnt++;
    }
    return p_cnt;
}

void init_parity(HAMMING_INFO_T *info)
{
    uint32_t n = info->n;
    uint32_t m = info->m;
    uint8_t **parity = info->parity;

    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t j = 0; j < m; j++)
        {
            parity[i][j] = 0;
        }
    }

    for (uint32_t i = 0; i < n; i++)
    {
        uint32_t shift = (1 << i);
        for (uint32_t j = shift - 1; j < m; j += (shift << 1))
        {
            uint32_t to = MIN(j + shift, m);
            for (uint32_t k = j; k < to; k++)
            {
                parity[i][k] = 1;
            }
        }
    }
}

void init_hash_table(HAMMING_INFO_T *info)
{
    uint32_t syndrome;
    for (uint32_t i = 0; i < (1 << info->n); i++)
    {
        info->hash_table[i] = UINT32_MAX;
    }
    for (uint32_t i = 0; i < info->m; i++)
    {
        syndrome = 0;
        for (uint32_t j = 0; j < info->n; j++)
        {
            syndrome <<= 1;
            syndrome |= info->parity[j][i];
        }
        assert(syndrome < (1 << info->n));
        info->hash_table[syndrome] = i;
    }
}

HAMMING_INFO_T *init_hamming(uint32_t len, uint32_t ded_on)
{
    HAMMING_INFO_T *info = (HAMMING_INFO_T *) malloc(sizeof(HAMMING_INFO_T));
    uint32_t p_cnt = cal_parity_count(len);
    info->ded_on = ded_on & 1;
    info->n = p_cnt;
    info->m = len + info->n;
    info->code_len = info->m + info->ded_on;
    info->hash_table = malloc((1 << info->n) * sizeof(uint32_t));
    info->parity = malloc(info->n * sizeof(uint8_t *));
    for (uint32_t i = 0; i < info->n; i++)
    {
        info->parity[i] = malloc(info->m * sizeof(uint8_t));
    }
    init_parity(info);
    init_hash_table(info);
//    show_matrix(info);
//    show_table(info);
    return info;
}

void hamming_codeword(HAMMING_INFO_T *info, const uint8_t *data, uint8_t *codeword)
{
    uint32_t mask = (1 << info->n) - 1;
    uint32_t len = info->m - info->n;
    uint32_t index = 1;
    uint32_t p = 0;
    uint32_t i = 0;
    while (i < len)
    {
        if (is_power_of_2(index))
        {
            index++;
            continue;
        }
        if (data[i] & 1)
        {
            p = p ^ (index & mask);
        }
        codeword[index - 1] = data[i];
        index++;
        i++;
    }

    index = 1;
    uint32_t l = info->n;
    for (i = 0; i < l; i++)
    {
        codeword[index - 1] = p & 1;
        index <<= 1;
        p >>= 1;
    }

    if (info->ded_on)
    {
        p = 0;
        for (i = 0; i < info->m; i++)
        {
            p ^= codeword[i];
        }
        codeword[info->m] = p;
    }
}

ERROR_INFO_T hamming_get_error_info(HAMMING_INFO_T *info, const uint8_t *code_word)
{
    ERROR_INFO_T err_info = {UNKNOWN_ERROR, 0, 0};
    uint32_t syndrome = 0;
    uint32_t v;
    for (uint32_t i = 0; i < info->n; i++)
    {
        v = 0;
        for (uint32_t j = 0; j < info->m; j++)
        {
            v ^= (info->parity[i][j] & code_word[j]);
        }
        syndrome <<= 1;
        syndrome |= v;
    }

    if (info->ded_on)
    {
        uint32_t last = 0;

        for (uint32_t i = 0; i < info->code_len; i++)
        {
            last ^= code_word[i];
        }

        if (last)
        {
            err_info.error_type = ONE_BIT_ERROR;
            err_info.error_index = syndrome == 0 ? info->m : info->hash_table[syndrome];
            err_info.error_syndrome = syndrome;
        }
        else
        {
            err_info.error_type = syndrome == 0 ? NO_BIT_ERROR : TWO_BIT_ERROR;
        }
    }
    else
    {
        err_info.error_type = syndrome == 0 ? NO_BIT_ERROR : ONE_BIT_ERROR;
        err_info.error_index = info->hash_table[syndrome];
        err_info.error_syndrome = syndrome;
    }

    return err_info;
}

uint32_t get_valid_error_count(HAMMING_INFO_T *info)
{
    return gen_random_number() % (2 + info->ded_on);
}

void inject_error(HAMMING_INFO_T *info, uint8_t *code_word, uint32_t error_count)
{
    if (error_count == 0 || error_count > 2 || (error_count == 2 && !info->ded_on))
    {
        return;
    }

    uint32_t index[2];
    index[0] = gen_random_number() % info->code_len;
    info->err_index = index[0];
    if (error_count == 2)
    {
        index[1] = gen_random_number() % info->code_len;
        while (index[1] == index[0])
        {
            index[1] = gen_random_number() % info->code_len;
        }
    }
    for (uint32_t i = 0; i < error_count; i++)
    {
        code_word[index[i]] ^= 1;
    }
}

int main()
{
    srand(time(NULL));
    uint32_t data_len = 256;
    uint32_t ded_on = 1;

    HAMMING_INFO_T *info = init_hamming(data_len, ded_on);
    uint8_t *data = malloc(data_len * sizeof(uint8_t));
    uint8_t *codeword = malloc(info->code_len * sizeof(uint8_t));
    uint32_t cnt[3] = {0, 0, 0};

    printf("n %d, m %d, code_len %d\n", info->n, info->m, info->code_len);

    for (uint32_t i = 0; i < 3000000; i++)
    {
        uint32_t num = gen_random_number() & ((1 << data_len) - 1);
        for (uint32_t j = 0; j < data_len; j++)
        {
            data[j] = num & 1;
            num >>= 1;
        }

        hamming_codeword(info, data, codeword);
        uint32_t error_count = get_valid_error_count(info);
        inject_error(info, codeword, error_count);
        ERROR_INFO_T err_info = hamming_get_error_info(info, codeword);

        switch (error_count)
        {
            case 0:
                assert(err_info.error_type == NO_BIT_ERROR);
                cnt[0]++;
                break;
            case 1:
                assert(err_info.error_type == ONE_BIT_ERROR);
                assert(err_info.error_index == info->err_index);
                cnt[1]++;
                break;
            case 2:
                assert(err_info.error_type == TWO_BIT_ERROR);
                cnt[2]++;
                break;
            default:
                break;
        }
    }

    printf("%d, %d, %d\n", cnt[0], cnt[1], cnt[2]);

    return 0;
}

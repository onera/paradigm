#include<stdio.h>
#include<stdlib.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_cuda_print.h"


int main(void) {
    printf("Hello World from host!\n");
    CUDA_print();
    return 0;
}

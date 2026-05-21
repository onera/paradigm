cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.16)

project(fake_install LANGUAGES C)

find_package(pdm
             REQUIRED COMPONENTS pdm pdmf
             HINTS $ENV{PDM_ROOT}
             )
/*
 * Copyright 2016 Thiago Nascimento nascimenthiago@gmail.com
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef __TEST_SUITE_MATRIX_H__
#define __TEST_SUITE_MATRIX_H__

#include <assert.h>
#include "../CommonFiles/matrix_parallel.h"

void run_all_test_matrix();
void test_parallel_max_wavefront();
void test_parallel_rms_wavefront();
void test_parallel_max_wavefront_apothen();

#endif /* __TEST_SUITE_MATRIX_H__ */

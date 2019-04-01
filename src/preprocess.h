/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford University. All Rights Reserved.
* Licensed under the “Non-exclusive Research Use License for SPRING Software” license (the "License");
* You may not use this file except in compliance with the License.
* The License is included in the distribution as license.pdf file.
 
* Software distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and limitations under the License.
*/

#ifndef SPRING_PREPROCESS_H_
#define SPRING_PREPROCESS_H_

#include <string>
#include "util.h"

namespace spring {

void preprocess(const std::string &infile_1, const std::string &infile_2,
                const std::string &temp_dir, compression_params &cp, const bool &gzip_flag);

}  // namespace spring

#endif  // SPRING_PREPROCESS_H_

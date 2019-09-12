#include <stdexcept>
#include <string>
#include "encoder.h"
#include "reorder.h"
#include "util.h"

namespace spring {

void call_reorder(const std::string &temp_dir, compression_params &cp) {
  size_t bitset_size_reorder = (2 * cp.max_readlen - 1) / 64 * 64 + 64;
  switch (bitset_size_reorder) {
    case 64:
      reorder_main<64>(temp_dir, cp);
      break;
    case 128:
      reorder_main<128>(temp_dir, cp);
      break;
    case 192:
      reorder_main<192>(temp_dir, cp);
      break;
    case 256:
      reorder_main<256>(temp_dir, cp);
      break;
    case 320:
      reorder_main<320>(temp_dir, cp);
      break;
    case 384:
      reorder_main<384>(temp_dir, cp);
      break;
    case 448:
      reorder_main<448>(temp_dir, cp);
      break;
    case 512:
      reorder_main<512>(temp_dir, cp);
      break;
    case 576:
      reorder_main<576>(temp_dir, cp);
      break;
    case 640:
      reorder_main<640>(temp_dir, cp);
      break;
    case 704:
      reorder_main<704>(temp_dir, cp);
      break;
    case 768:
      reorder_main<768>(temp_dir, cp);
      break;
    case 832:
      reorder_main<832>(temp_dir, cp);
      break;
    case 896:
      reorder_main<896>(temp_dir, cp);
      break;
    case 960:
      reorder_main<960>(temp_dir, cp);
      break;
    case 1024:
      reorder_main<1024>(temp_dir, cp);
      break;
    default:
      throw std::runtime_error("Wrong bitset size.");
  }
}

void call_encoder(const std::string &temp_dir, compression_params &cp) {
  size_t bitset_size_encoder = (3 * cp.max_readlen - 1) / 64 * 64 + 64;
  switch (bitset_size_encoder) {
    case 64:
      encoder_main<64>(temp_dir, cp);
      break;
    case 128:
      encoder_main<128>(temp_dir, cp);
      break;
    case 192:
      encoder_main<192>(temp_dir, cp);
      break;
    case 256:
      encoder_main<256>(temp_dir, cp);
      break;
    case 320:
      encoder_main<320>(temp_dir, cp);
      break;
    case 384:
      encoder_main<384>(temp_dir, cp);
      break;
    case 448:
      encoder_main<448>(temp_dir, cp);
      break;
    case 512:
      encoder_main<512>(temp_dir, cp);
      break;
    case 576:
      encoder_main<576>(temp_dir, cp);
      break;
    case 640:
      encoder_main<640>(temp_dir, cp);
      break;
    case 704:
      encoder_main<704>(temp_dir, cp);
      break;
    case 768:
      encoder_main<768>(temp_dir, cp);
      break;
    case 832:
      encoder_main<832>(temp_dir, cp);
      break;
    case 896:
      encoder_main<896>(temp_dir, cp);
      break;
    case 960:
      encoder_main<960>(temp_dir, cp);
      break;
    case 1024:
      encoder_main<1024>(temp_dir, cp);
      break;
    case 1088:
      encoder_main<1088>(temp_dir, cp);
      break;
    case 1152:
      encoder_main<1152>(temp_dir, cp);
      break;
    case 1216:
      encoder_main<1216>(temp_dir, cp);
      break;
    case 1280:
      encoder_main<1280>(temp_dir, cp);
      break;
    case 1344:
      encoder_main<1344>(temp_dir, cp);
      break;
    case 1408:
      encoder_main<1408>(temp_dir, cp);
      break;
    case 1472:
      encoder_main<1472>(temp_dir, cp);
      break;
    case 1536:
      encoder_main<1536>(temp_dir, cp);
      break;
    default:
      throw std::runtime_error("Wrong bitset size.");
  }
}
}  // namespace spring

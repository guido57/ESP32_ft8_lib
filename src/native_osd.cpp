#include "native_osd.h"

#include "ft8/constants.h"
#include "ft8/crc.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

extern "C" void encode174(const uint8_t* message, uint8_t* codeword);

namespace {

constexpr int kBits = FTX_LDPC_N;
constexpr int kMessageBits = FTX_LDPC_K;
constexpr int kWeakPositions = 10;
#ifndef FT8_NATIVE_OSD_MAX_FLIPS
#define FT8_NATIVE_OSD_MAX_FLIPS 2
#endif
constexpr int kMaxFlips = FT8_NATIVE_OSD_MAX_FLIPS;
using Row = std::bitset<kMessageBits>;

bool packed_bit(const uint8_t bytes[], int bit)
{
    return (bytes[bit / 8] & (uint8_t(0x80) >> (bit % 8))) != 0;
}

void set_packed_bit(uint8_t bytes[], int bit)
{
    bytes[bit / 8] |= uint8_t(0x80) >> (bit % 8);
}

bool crc_ok(const uint8_t message[])
{
    uint8_t input[FTX_LDPC_K_BYTES];
    std::memcpy(input, message, sizeof(input));
    const uint16_t expected = ftx_extract_crc(input);
    input[9] &= 0xf8u;
    input[10] = 0;
    return expected == ftx_compute_crc(input, 82);
}

void codeword_to_bits(const uint8_t codeword[], uint8_t bits[])
{
    for (int i = 0; i < kBits; ++i)
        bits[i] = packed_bit(codeword, i) ? 1 : 0;
}

double metric(const float log174[], const uint8_t codeword[])
{
    double sum = 0.0;
    for (int bit = 0; bit < kBits; ++bit)
        sum += packed_bit(codeword, bit) ? log174[bit] : -log174[bit];
    return sum;
}

bool invert_binary_matrix(const std::array<Row, kMessageBits>& input,
                          std::array<Row, kMessageBits>* inverse)
{
    std::array<std::bitset<2 * kMessageBits>, kMessageBits> augmented{};
    for (int row = 0; row < kMessageBits; ++row) {
        for (int col = 0; col < kMessageBits; ++col)
            augmented[row][col] = input[row][col];
        augmented[row][kMessageBits + row] = true;
    }

    for (int col = 0; col < kMessageBits; ++col) {
        int pivot = col;
        while (pivot < kMessageBits && !augmented[pivot][col])
            ++pivot;
        if (pivot == kMessageBits)
            return false;
        std::swap(augmented[col], augmented[pivot]);
        for (int row = 0; row < kMessageBits; ++row) {
            if (row != col && augmented[row][col])
                augmented[row] ^= augmented[col];
        }
    }

    for (int row = 0; row < kMessageBits; ++row)
        for (int col = 0; col < kMessageBits; ++col)
            (*inverse)[row][col] = augmented[row][kMessageBits + col];
    return true;
}

void solve_message(const std::array<Row, kMessageBits>& inverse,
                   const std::array<uint8_t, kMessageBits>& y,
                   uint8_t message[])
{
    std::memset(message, 0, FTX_LDPC_K_BYTES);
    for (int row = 0; row < kMessageBits; ++row) {
        bool value = false;
        for (int col = 0; col < kMessageBits; ++col)
            value ^= inverse[row][col] && y[col];
        if (value)
            set_packed_bit(message, row);
    }
}

} // namespace

extern "C" int native_osd_recover(const float log174[], uint8_t recovered174[], float* score)
{
    // Derive the systematic generator once from the project's encoder.  This
    // avoids embedding a second copy of the FT8 parity matrix in the native
    // decoder and makes the list decoder follow the same on-air code.
    std::array<Row, kBits> generator{};
    for (int message_bit = 0; message_bit < kMessageBits; ++message_bit) {
        uint8_t message[FTX_LDPC_K_BYTES]{};
        uint8_t codeword[FTX_LDPC_N_BYTES]{};
        set_packed_bit(message, message_bit);
        encode174(message, codeword);
        for (int code_bit = 0; code_bit < kBits; ++code_bit)
            generator[code_bit][message_bit] = packed_bit(codeword, code_bit);
    }

    std::vector<int> order(kBits);
    for (int i = 0; i < kBits; ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [log174](int a, int b) {
        return std::fabs(log174[a]) > std::fabs(log174[b]);
    });

    // Pick the most reliable linearly independent observations.
    std::array<Row, kMessageBits> echelon{};
    std::array<Row, kMessageBits> selected_rows{};
    std::array<int, kMessageBits> selected_bits{};
    int rank = 0;
    for (int code_bit : order) {
        Row trial = generator[code_bit];
        for (int pivot = 0; pivot < kMessageBits; ++pivot) {
            if (!trial[pivot])
                continue;
            if (echelon[pivot].none()) {
                echelon[pivot] = trial;
                selected_rows[rank] = generator[code_bit];
                selected_bits[rank] = code_bit;
                ++rank;
                break;
            }
            trial ^= echelon[pivot];
        }
        if (rank == kMessageBits)
            break;
    }
    if (rank != kMessageBits)
        return 0;

    std::array<Row, kMessageBits> inverse{};
    if (!invert_binary_matrix(selected_rows, &inverse))
        return 0;

    std::array<uint8_t, kMessageBits> y{};
    for (int i = 0; i < kMessageBits; ++i)
        y[i] = log174[selected_bits[i]] > 0.0f;

    std::array<int, kMessageBits> weak_order{};
    for (int i = 0; i < kMessageBits; ++i)
        weak_order[i] = i;
    std::sort(weak_order.begin(), weak_order.end(), [&](int a, int b) {
        return std::fabs(log174[selected_bits[a]]) < std::fabs(log174[selected_bits[b]]);
    });

    bool found = false;
    double best_metric = 0.0;
    uint8_t best_codeword[FTX_LDPC_N_BYTES]{};
    const int weak_count = std::min(kWeakPositions, kMessageBits);
    std::array<int, kMaxFlips> flips{};
    auto evaluate = [&](int flip_count) {
        auto trial_y = y;
        for (int i = 0; i < flip_count; ++i)
            trial_y[weak_order[flips[i]]] ^= 1;

        uint8_t message[FTX_LDPC_K_BYTES]{};
        solve_message(inverse, trial_y, message);
        if (!crc_ok(message))
            return;

        uint8_t codeword[FTX_LDPC_N_BYTES]{};
        encode174(message, codeword);
        const double candidate_metric = metric(log174, codeword);
        if (!found || candidate_metric > best_metric) {
            std::memcpy(best_codeword, codeword, sizeof(best_codeword));
            best_metric = candidate_metric;
            found = true;
        }
    };
    auto enumerate = [&](auto&& self, int first, int depth, int flip_count) -> void {
        if (depth == flip_count) {
            evaluate(flip_count);
            return;
        }
        for (int i = first; i <= weak_count - (flip_count - depth); ++i) {
            flips[depth] = i;
            self(self, i + 1, depth + 1, flip_count);
        }
    };
    for (int flip_count = 0; flip_count <= kMaxFlips; ++flip_count)
        enumerate(enumerate, 0, 0, flip_count);

    if (!found)
        return 0;
    codeword_to_bits(best_codeword, recovered174);
    if (score != nullptr)
        *score = static_cast<float>(best_metric);
    return 1;
}

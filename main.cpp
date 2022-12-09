#include <iostream>
#include <vector>
#include <numeric>
#include <queue>
#include <map>
#include <functional>
#include <set>
#include <thread>
#include <future>
#include <cmath>
#include <cassert>

#define NUM_OF_THREADS 4

/*
 The code in this file finds the master key for the wes/wes.c file where the master key is hardcoded.

 Compile/Run instructions:
    g++ -O3 main.cpp -o wes-cryptanalysis
    ./wes-cryptanalysis
 */

// S-box Table
int sBox[8][16] = {
        /* 1 */
        // 0, 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  13, 14, 15
        { 6,  12, 3,  8,  14, 5,  11, 1,  2,  4,  13, 7,  0,  10, 15, 9 },
        /* 2 */
        { 10, 14, 15, 11, 6,  8,  3,  13, 7,  9,  2,  12, 1,  0,  4,  5 },
        /* 3 */
        { 10, 0,  9,  14, 6,  3,  15, 5,  1,  13, 12, 7,  11, 4,  2,  8 },
        /* 4 */
        { 15, 9,  7,  0,  10, 13, 2,  4,  3,  6,  12, 5,  1,  8,  14, 11 },
        /* 5 */
        { 2,  12, 4,  1,  7,  10, 11, 6,  8,  5,  3,  15, 13, 0,  14, 9 },
        /* 6 */
        { 12, 1,  10, 15, 9,  2,  6,  8,  0,  13, 3,  4,  14, 7,  5,  11 },
        /* 7 */
        { 4,  11, 2,  14, 15, 0,  8,  13, 3,  12, 9,  7,  5,  10, 6,  1 },
        /* 8 */
        { 13, 2,  8,  4,  6,  15, 11, 1,  10, 9,  3,  14, 5,  0,  12, 7 }
};

// Permutation Table
int per[32] = {
        16, 7, 20, 21,
        29, 12, 28, 17,
        1, 15, 23, 26,
        5, 18, 31, 10,
        2, 8, 24, 14,
        32, 27, 3, 9,
        19, 13, 30, 6,
        22, 11, 4, 25
};

uint32_t mask[32]; /* permutation mask to speed up the permutation transform */

uint32_t sbox_layer( uint32_t x )
{
    uint32_t res = 0;
    res = res | ( sBox[ 0 ][ ( x >> 28 ) & 0xf ] << 28 );
    res = res | ( sBox[ 1 ][ ( x >> 24 ) & 0xf ] << 24 );
    res = res | ( sBox[ 2 ][ ( x >> 20 ) & 0xf ] << 20 );
    res = res | ( sBox[ 3 ][ ( x >> 16 ) & 0xf ] << 16 );
    res = res | ( sBox[ 4 ][ ( x >> 12 ) & 0xf ] << 12 );
    res = res | ( sBox[ 5 ][ ( x >> 8 ) & 0xf ] << 8 );
    res = res | ( sBox[ 6 ][ ( x >> 4 ) & 0xf ] << 4 );
    res = res | ( sBox[ 7 ][ x & 0xf ] );
    return res;
}

uint32_t permute( uint32_t x )
{
    uint32_t res = 0;
    for ( int i = 0; i < 32; i++ )
        res |= ( ( x & mask[ i ] ) << ( per[ i ] - 1 ) ) >> i;
    return res;
}

/**
 *  WES round function:
 *
 *          1) xor with the round key
 *          2) Pass through f-box:
 *             -- s-boxes
 *             -- permutation
 *
 *
 *                         +------------- K (in)
 *                         |
 *           +------+      v
 *    out <--|  f   |<--- xor <--- x (in)
 *           +------+
 *
 *
 *  f(x) : out <-- PERMUTATION_BOX <-- SBOX's <-- x
 *
 * */
uint32_t round_func( uint32_t x, uint32_t rkey )
{
    x = x ^ rkey;
    x = sbox_layer( x );
    x = permute( x );
    return x;
}

/* Optimization: mask is used to extract bits at certain position.
 * Since we can do it once before the encryption, it will save us
 * some operations during the encryption */
int precompute_wes_permutation_mask()
{
    for ( int i = 0; i < 32; i++ )
        mask[ i ] = 1 << ( 32 - per[ i ] );
    return 0;

}

void print2DMatrix( const std::vector< std::vector< int>>& differentialTable )
{
    std::cout << "[\n";
    for ( int i = 0; i < differentialTable.size(); ++i )
    {
        for ( int j = 0; j < differentialTable[ i ].size(); ++j )
        {
            std::cout << differentialTable[ i ][ j ] << ", ";
        }
        std::cout << "\n";
    }
    std::cout << "]\n";
}

/* Function that computes the difference distribution tables for s-boxes 1 to 8 */
std::vector< std::vector< std::vector< int>> > sBoxDifferenceDistributionTables()
{
    std::vector< std::vector< std::vector< int>> > sBoxTables( 8 );
    for ( int sBoxNbr = 0; sBoxNbr < 8; ++sBoxNbr )
    {
        std::vector< std::vector< int>> differentialTable( 16, std::vector< int >( 16, 0 ) );

        for ( int deltaX = 0; deltaX < 16; ++deltaX )
        {
            for ( int x1 = 0; x1 < 16; ++x1 )
            {
                int x2 = x1 ^ deltaX;
                int y2 = sBox[ sBoxNbr ][ x2 ];
                int y1 = sBox[ sBoxNbr ][ x1 ];
                int deltaY = y1 ^ y2;
                ++differentialTable[ deltaX ][ deltaY ]; // increment the occurrence of the output differential given the input differential
            }
        }

        std::cout << "Differential table for s-box: " << sBoxNbr + 1 << std::endl;
        print2DMatrix( differentialTable );
        sBoxTables[ sBoxNbr ] = differentialTable;
    }
    return sBoxTables;
}

std::vector< std::vector< std::vector< std::pair< uint32_t, int > >> >
optimizeDifferenceDistributionTables( const std::vector< std::vector< std::vector< int>> >& diffTables )
{
    std::vector< std::vector< std::vector< std::pair< uint32_t, int > >> > optimizedDiffTables( diffTables.size() );

    for ( int sBoxNbr = 0; sBoxNbr < 8; ++sBoxNbr )
    {
        std::vector< std::vector< std::pair< uint32_t, int > >> matrix( 16 );
        for ( int deltaX = 0; deltaX < 16; ++deltaX )
        {
            std::vector< std::pair< uint32_t, int > > deltaXRow;
            for ( int deltaY = 0; deltaY < 16; ++deltaY )
            {
                if ( diffTables[ sBoxNbr ][ deltaX ][ deltaY ] != 0 )
                {
                    deltaXRow.emplace_back( std::pair( deltaY, diffTables[ sBoxNbr ][ deltaX ][ deltaY ] ) );
                }
            }
            matrix[ deltaX ] = deltaXRow;
        }
        optimizedDiffTables[ sBoxNbr ] = matrix;
    }

    return optimizedDiffTables;
}

// Function that each thread will execute to break subkey k4
std::set< uint32_t > workerThreadToCrackKey4( uint32_t begin,
                                              uint32_t end,
                                              std::vector< std::vector< uint64_t > > cipherPairs,
                                              uint32_t lastRoundInputDifference
)
{
    std::set< uint32_t > set;
    for ( uint32_t k4 = begin; k4 < end; ++k4 )
    {
        int score = 0;

        for ( const std::vector< uint64_t >& pairs : cipherPairs )
        {
            uint64_t c1 = pairs[ 1 ];
            uint64_t c2 = pairs[ 3 ];
            uint64_t deltaC = ( c1 ^ c2 );

            uint32_t leftHalf = deltaC >> 32;
            uint32_t rightHalfC1 = c1 & 0xffffffff;
            uint32_t rightHalfC2 = c2 & 0xffffffff;

            uint32_t result1 = round_func( rightHalfC1, k4 );
            uint32_t result2 = round_func( rightHalfC2, k4 );
            uint32_t outputDifference = result1 ^ result2;

            uint32_t calculatedInputDifference = outputDifference ^ leftHalf;

            if ( calculatedInputDifference == lastRoundInputDifference )
            {
                ++score;
            }
        }

        if ( score == cipherPairs.size() )
        {
            // std::cout << "Found subkey k4: " << " 0x" << k4 << std::endl;
            set.insert( k4 );
        }
    }

    std::cout << std::dec;
    std::cout << "There are " << set.size() << " matching subkeys\n";
    return set;
}

std::set< uint32_t > crackKey4()
{
    // Let's try to break the last round partial key k4

    // [Plaintext 1, Ciphertext 1, P2, C2]
    std::vector< std::vector< uint64_t > > deltaPairs = {
            { 0xe0a621281e76a643, 0xFF3CCAC4EE3FFC8D, 0x82a021281e76a643, 0x264CD21D883A7D2C },
            { 0x9d5da5a9eff0f367, 0xF735CA49CF34C3BD, 0xff5ba5a9eff0f367, 0x4AEA5CBF89A04C9F },
            { 0x1629032bf2d0f365, 0x19CCF3B531745F36, 0x742f032bf2d0f365, 0xBBF9A2B05FE3D496 },
            { 0xf2265d81fc1bd964, 0xBE9D1520BEA28942, 0x90205d81fc1bd964, 0x61E828D0D8260AE3 },
            { 0x03a54c0a20e8f05d, 0x59813499F0F5B72C, 0x61a34c0a20e8f05d, 0xCA89424AB5F5360C },
            { 0x237f61071881097b, 0xCA2FD612D8659769, 0x417961071881097b, 0x9E36963BB4E39EEC },
            { 0x44a3c3ecbdb612a2, 0xD6F435EF0712C355, 0x26a5c3ecbdb612a2, 0xE1CBDD10429348F7 },
            { 0x1f27d48d2c0fe9c8, 0x6EF64F5918F4138C, 0x7d21d48d2c0fe9c8, 0x88CB5BCD76E39E0C }
    };

    // We know that the input difference to the last round
    uint32_t lastRoundInputDifference = 0x8080D052;
    std::cout << std::hex; // numbers will be output in hex format

    std::vector< std::packaged_task< std::set< uint32_t >( uint32_t,
                                                           uint32_t,
                                                           std::vector< std::vector< uint64_t > >,
                                                           uint32_t
    ) > > tasksForThreads( NUM_OF_THREADS );

    std::vector< std::future< std::set< uint32_t > > > futures;
    std::vector< std::thread > threads;

    int increment = ( int ) std::pow( 2, 32 ) / NUM_OF_THREADS;
    int begin = 0;
    int end = increment;

    for ( int i = 0; i < NUM_OF_THREADS; ++i )
    {
        std::packaged_task< std::set< uint32_t >( uint32_t,
                                                  uint32_t,
                                                  std::vector< std::vector< uint64_t > >,
                                                  uint32_t
        ) > task { workerThreadToCrackKey4 }; // Create task
        futures.emplace_back( task.get_future() );      // Get the future object.

        // std::packaged_task is move-only
        threads.emplace_back( std::move( task ), begin, end, deltaPairs, lastRoundInputDifference );

        begin = end;
        end += increment;
    }

    std::cout << "A total of " << NUM_OF_THREADS << " threads have started to break round 4...\n";

    // wait for all threads to finish
    for ( std::thread& t : threads )
    {
        t.join();
    }

    std::cout << "All threads have finished their work for breaking sub-key k4...\n";

    // Extract the result of each thread
    std::set< uint32_t > subKeysK4;

    for ( int i = 0; i < futures.size(); i += 2 )
    {
        std::set< uint32_t > set1 = futures[ i ].get();
        std::set< uint32_t > set2 = futures[ i + 1 ].get();

        std::set_union(
                set1.begin(), set1.end(),
                set2.begin(), set2.end(),
                std::inserter( subKeysK4, std::begin( subKeysK4 ) )
        );
    }

    std::cout << "Final set has " << subKeysK4.size() << " candidate sub-keys k4\n";
    return subKeysK4;
}

struct candidateKeysK3K4
{
    std::set< uint32_t > candidateK3s;
    std::set< uint32_t > candidateK4s;
};

candidateKeysK3K4 workerThreadToCrackKey3(
        uint32_t begin,
        uint32_t end,
        std::set< uint32_t > subKeysK4,
        std::vector< std::vector< uint64_t > > pairs,
        uint32_t diffCharacteristicToRound3
)
{
    std::set< uint32_t > candidateK3s;
    std::set< uint32_t > candidateK4s;

    for ( uint32_t leftMostBitsOfK3 = begin; leftMostBitsOfK3 < end; ++leftMostBitsOfK3 )
    {
        for ( uint32_t k4 : subKeysK4 )
        {
            int score = 0;

            // Half of key k3 is already present in our candidate key k4, we can simply extract those bits
            uint16_t rightMostBitsOfK3 = 0x0;

            int bitPosition = 30;
            int shiftAmount = 15;
            for ( int i = 0; i < 16; ++i )
            {
                // Take the bit that is shared between k3 and k4
                uint16_t correspondingBitInK4 = ( k4 & ( 1 << bitPosition ) ) > 0 ? 1 : 0;
                rightMostBitsOfK3 |= ( correspondingBitInK4 << shiftAmount );
                bitPosition -= 2;
                --shiftAmount;
            }

            uint32_t k3 = ( leftMostBitsOfK3 << 16 ) | rightMostBitsOfK3;

            for ( const std::vector< uint64_t >& pair : pairs )
            {
                uint64_t c1 = pair[ 1 ];
                uint64_t c2 = pair[ 3 ];

                uint32_t rightHalfC1 = c1 & 0xffffffff;
                uint32_t rightHalfC2 = c2 & 0xffffffff;
                uint32_t leftHalfC1 = c1 >> 32;
                uint32_t leftHalfC2 = c2 >> 32;

                uint64_t deltaC = ( c1 ^ c2 );
                uint32_t rightHalf = deltaC & 0xffffffff;

                uint32_t result1 = round_func( rightHalfC1, k4 ) ^ leftHalfC1;
                uint32_t result2 = round_func( rightHalfC2, k4 ) ^ leftHalfC2;

                uint32_t round3Result1 = round_func( result1, k3 );
                uint32_t round3Result2 = round_func( result2, k3 );

                uint32_t out = round3Result1 ^ round3Result2 ^ rightHalf;

                if ( diffCharacteristicToRound3 == out )
                {
                    ++score;
                }
            }

            if ( score == pairs.size() )
            {
                candidateK3s.insert( k3 );
                candidateK4s.insert( k4 );
            }
        }
    }
    return candidateKeysK3K4 { candidateK3s, candidateK4s };
}

std::pair< std::set< uint32_t >, std::set< uint32_t > > crackKey3( const std::set< uint32_t >& candidateKeysK4 )
{
    // Let's break sub-key k3
    std::cout << "\nFinished loading sub-keys into set...proceeding with cracking round 3\n";

    // [Plaintext 1, Ciphertext 1, P2, C2]
    std::vector< std::vector< uint64_t > > pairs = {
            { 0x4a28d313c14be57e, 0x0A1EA4A827794303, 0xb5d72cecc14be57e, 0x84EE886538EA684D },
            { 0xfb870040f62f8a43, 0x12DF803607E99015, 0x0478ffbff62f8a43, 0xACBC31D5C88550BD },
            { 0xfaa4aa317d1f419c, 0x5D527DA3F87C07F5, 0x055b55ce7d1f419c, 0xFD389A200E034942 },
            { 0xbfdd81549096fdd5, 0x1B39BCBF5FFF6F7C, 0x40227eab9096fdd5, 0xBAEDE607106F5B4E },
            { 0x8d54144754bffc8b, 0xBD240FFAE2C28D7B, 0x72abebb854bffc8b, 0xAB7E6D6B4AC39E18 },
            { 0xa6f98171d7d9e1ff, 0xC0FA679AF96E3124, 0x59067e8ed7d9e1ff, 0xA1E77E206ABE1926 },
            { 0x49b6478779b82e5b, 0xD53599464D1A5D5E, 0xb649b87879b82e5b, 0x368D0E6E3FE1EA42 },
            { 0xd4184e97c7c6b12e, 0x680FB9A85CF87F02, 0x2be7b168c7c6b12e, 0x9C76031A70884878 },
            { 0xcaeaac6eb1e1e847, 0x3E7F118BCD53E0B7, 0x35155391b1e1e847, 0x68D7454A6354E961 },
            { 0xd3ef23f21c62d8fc, 0x8831A458417E9A64, 0x2c10dc0d1c62d8fc, 0x3272A8A9B0492ABB },
            { 0x25dea086b2c71305, 0xCD2BC98DA8256ED0, 0xda215f79b2c71305, 0x1743FE4C0AB147D3 },
            { 0x55b91ec107bef453, 0x709FBA067566D69B, 0xaa46e13e07bef453, 0x42888F2CFABC9A82 },
            { 0xa6c79ffea4b986f5, 0x02F8C59885C24040, 0x59386001a4b986f5, 0x6051FB3377C7F4D7 }
    };

    uint32_t diffCharacteristicToRound3 = 0xffffffff;
    std::cout << std::hex;

    std::vector< std::packaged_task< std::set< uint32_t >( uint32_t,
                                                           uint32_t,
                                                           std::vector< std::vector< uint64_t > >,
                                                           uint32_t
    ) > > tasksForThreads( NUM_OF_THREADS );

    std::vector< std::future< candidateKeysK3K4 > > futures;
    std::vector< std::thread > threads;

    int increment = ( int ) std::pow( 2, 16 ) / NUM_OF_THREADS;
    int begin = 0;
    int end = increment;

    for ( int i = 0; i < NUM_OF_THREADS; ++i )
    {
        std::packaged_task< candidateKeysK3K4( uint32_t begin,
                                               uint32_t end,
                                               std::set< uint32_t > subKeysK4,
                                               std::vector< std::vector< uint64_t > > pairs,
                                               uint32_t diffCharacteristicToRound3
        ) > task { workerThreadToCrackKey3 }; // Create task
        futures.emplace_back( task.get_future() );      // Get the future object.

        // std::packaged_task is move-only
        threads.emplace_back( std::move( task ), begin, end, candidateKeysK4, pairs, diffCharacteristicToRound3 );

        begin = end;
        end += increment;
    }

    std::cout << "A total of " << NUM_OF_THREADS << " threads have started to break round 3...\n";

    // wait for all threads to finish
    for ( std::thread& t : threads )
    {
        t.join();
    }

    std::cout << "All threads have finished their work for breaking sub-key k3...\n";

    // Extract the result of each thread
    std::set< uint32_t > k4;
    std::set< uint32_t > k3;

    for ( int i = 0; i < futures.size(); i += 2 )
    {
        candidateKeysK3K4 keys1 = futures[ i ].get();
        candidateKeysK3K4 keys2 = futures[ i + 1 ].get();

        std::set_union(
                keys1.candidateK3s.begin(), keys1.candidateK3s.end(),
                keys2.candidateK3s.begin(), keys2.candidateK3s.end(),
                std::inserter( k3, std::begin( k3 ) )
        );

        std::set_union(
                keys1.candidateK4s.begin(), keys1.candidateK4s.end(),
                keys2.candidateK4s.begin(), keys2.candidateK4s.end(),
                std::inserter( k4, std::begin( k4 ) )
        );
    }

    std::cout << "K4 candidate keys: \n";
    for ( uint32_t k : k4 )
    {
        std::cout << k << std::endl;
    }

    std::cout << "\nK3 candidate keys: \n";
    for ( uint32_t k : k3 )
    {
        std::cout << k << std::endl;
    }

    return { k3, k4 };
}

struct candidateKeysK2K3K4
{
    std::set< uint32_t > candidateK2s;
    std::set< uint32_t > candidateK3s;
    std::set< uint32_t > candidateK4s;
};

candidateKeysK2K3K4 workerThreadToCrackKey2(
        uint32_t begin,
        uint32_t end,
        std::set< uint32_t > subKeysK3,
        std::set< uint32_t > subKeysK4,
        std::vector< std::vector< uint64_t > > plaintextPairs
)
{
    std::set< uint32_t > candidateK2s;
    std::set< uint32_t > candidateK3s;
    std::set< uint32_t > candidateK4s;

    for ( uint32_t k2 = begin; k2 < end; ++k2 )
    {
        for ( uint32_t k3 : subKeysK3 )
        {
            // half of the bits of k2 are already present in k3, so we can simply check if the bits between k2 and k3 match
            uint16_t leftMostBitsOfK2 = 0x0;

            int bitPosition = 30;
            int shiftAmount = 15;
            for ( int i = 0; i < 16; ++i )
            {
                // Take the bit that is shared between k2 and k3
                uint16_t bitInK2 = ( k2 & ( 1 << bitPosition ) ) > 0 ? 1 : 0;
                leftMostBitsOfK2 |= ( bitInK2 << shiftAmount );
                bitPosition -= 2;
                --shiftAmount;
            }

            // First 16 bits of k2 that should be in k3
            uint16_t firstHalfOfK3 = k3 >> 16;

            if ( leftMostBitsOfK2 != firstHalfOfK3 ) continue; // shared bits between k2 and k3 don't match

            for ( uint32_t k4 : subKeysK4 )
            {
                int score = 0;
                for ( const std::vector< uint64_t >& plaintext : plaintextPairs )
                {
                    uint32_t L1 = plaintext[ 0 ] & 0xffffffff;
                    uint32_t CL = plaintext[ 1 ] >> 32;
                    uint32_t CR = plaintext[ 1 ] & 0xffffffff;

                    uint32_t R2 = round_func( CR, k4 ) ^ CL;
                    uint32_t L2 = round_func( R2, k3 ) ^ CR;

                    uint32_t R1 = L2;
                    uint32_t L1Backwards = round_func( R1, k2 ) ^ R2;

                    if ( L1Backwards == L1 )
                    {
                        ++score; // we have a match for k2
                    }
                }

                if ( score == plaintextPairs.size() )
                {
                    std::cout << "Found a matching k2! -> 0x" << k2 << std::endl;
                    candidateK2s.insert( k2 );
                    candidateK3s.insert( k3 );
                    candidateK4s.insert( k4 );
                }
            }
        }
    }

    return candidateKeysK2K3K4 { candidateK2s, candidateK3s, candidateK4s };
}

void crackKey2( const std::set< uint32_t >& subKeysK3, const std::set< uint32_t >& subKeysK4 )
{
    // [Plaintext 1, Ciphertext 1]
    std::vector< std::vector< uint64_t > > plaintexts = {
            { 0x4a28d313c14be57e, 0x0A1EA4A827794303 },
            { 0xfb870040f62f8a43, 0x12DF803607E99015 },
            { 0xfaa4aa317d1f419c, 0x5D527DA3F87C07F5 },
            { 0xbfdd81549096fdd5, 0x1B39BCBF5FFF6F7C },
            { 0x8d54144754bffc8b, 0xBD240FFAE2C28D7B },
            { 0xa6f98171d7d9e1ff, 0xC0FA679AF96E3124 },
            { 0x49b6478779b82e5b, 0xD53599464D1A5D5E },
            { 0xd4184e97c7c6b12e, 0x680FB9A85CF87F02 },
            { 0xcaeaac6eb1e1e847, 0x3E7F118BCD53E0B7 },
            { 0xd3ef23f21c62d8fc, 0x8831A458417E9A64 },
            { 0x25dea086b2c71305, 0xCD2BC98DA8256ED0 },
            { 0x55b91ec107bef453, 0x709FBA067566D69B },
            { 0xa6c79ffea4b986f5, 0x02F8C59885C24040 }
    };

    std::cout << std::hex;

    std::vector< std::packaged_task< std::set< uint32_t >( uint32_t,
                                                           uint32_t,
                                                           std::vector< std::vector< uint64_t > >
    ) > > tasksForThreads( NUM_OF_THREADS );

    std::vector< std::future< candidateKeysK2K3K4 > > futures;
    std::vector< std::thread > threads;

    int increment = ( int ) std::pow( 2, 32 ) / NUM_OF_THREADS;
    int begin = 0;
    int end = increment;

    for ( int i = 0; i < NUM_OF_THREADS; ++i )
    {
        std::packaged_task< candidateKeysK2K3K4(
                uint32_t,
                uint32_t,
                std::set< uint32_t >,
                std::set< uint32_t >,
                std::vector< std::vector< uint64_t > >
        ) > task { workerThreadToCrackKey2 }; // Create task

        futures.emplace_back( task.get_future() );      // Get the future object.

        // std::packaged_task is move-only
        threads.emplace_back( std::move( task ), begin, end, subKeysK3, subKeysK4, plaintexts );

        begin = end;
        end += increment;
    }

    std::cout << "A total of " << NUM_OF_THREADS << " threads have started to break round 2...\n";

    // wait for all threads to finish
    for ( std::thread& t : threads )
    {
        t.join();
    }

    std::cout << "\nAll threads have finished their work for breaking sub-key k2...\n";

    // Extract the result of each thread
    std::set< uint32_t > k4;
    std::set< uint32_t > k3;
    std::set< uint32_t > k2;

    for ( int i = 0; i < futures.size(); i += 2 )
    {
        candidateKeysK2K3K4 keys1 = futures[ i ].get();
        candidateKeysK2K3K4 keys2 = futures[ i + 1 ].get();

        std::set_union(
                keys1.candidateK2s.begin(), keys1.candidateK2s.end(),
                keys2.candidateK2s.begin(), keys2.candidateK2s.end(),
                std::inserter( k2, std::begin( k2 ) )
        );

        std::set_union(
                keys1.candidateK3s.begin(), keys1.candidateK3s.end(),
                keys2.candidateK3s.begin(), keys2.candidateK3s.end(),
                std::inserter( k3, std::begin( k3 ) )
        );

        std::set_union(
                keys1.candidateK4s.begin(), keys1.candidateK4s.end(),
                keys2.candidateK4s.begin(), keys2.candidateK4s.end(),
                std::inserter( k4, std::begin( k4 ) )
        );
    }

    std::cout << "\nK2 candidate keys: \n";
    for ( uint32_t k : k2 )
    {
        std::cout << k << std::endl;
    }

    std::cout << "\nK4 candidate keys: \n";
    for ( uint32_t k : k4 )
    {
        std::cout << k << std::endl;
    }

    assert( k2.size() == 1 );
    assert( k4.size() == 1 );

    uint64_t firstHalf = *k2.begin();
    uint64_t masterKey = ( firstHalf << 32 ) | *k4.begin();

    std::cout << std::hex << "Found master key! -> 0x" << masterKey << std::endl;
}

int main()
{
    // Part 1: Compute difference distribution tables for all s-boxes
    precompute_wes_permutation_mask(); /* Just an optimization: makes permutation step a bit faster */

    std::vector< std::vector< std::vector< int>> > diffTables = sBoxDifferenceDistributionTables();
    std::vector< std::vector< std::vector< std::pair< uint32_t, int > >> > optimizedDiffTables = optimizeDifferenceDistributionTables(
            diffTables );

    std::set< uint32_t > candidateKeysK4 = crackKey4();
    std::pair< std::set< uint32_t >, std::set< uint32_t > > k3AndK4 = crackKey3( candidateKeysK4 );
    crackKey2( k3AndK4.first, k3AndK4.second );

    return 0;
}

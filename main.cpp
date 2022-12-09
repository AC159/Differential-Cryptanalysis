#include <iostream>
#include <vector>
#include <numeric>
#include <queue>
#include <map>
#include <functional>
#include <set>
#include <thread>
#include <future>
#include <fstream>

#define FINAL_ROUND 4

struct highProbabilityDifferential
{
    double probability;
    std::vector< std::pair< uint32_t, uint32_t > > inputsAndOutputs;
};

auto lambda = []( const highProbabilityDifferential& diff1, const highProbabilityDifferential& diff2 ) -> bool
{
    if ( diff1.probability <= diff2.probability ) return true;
    return false;
};

// Priority queue that contains all probabilities. The first element will be the highest probability
std::priority_queue< highProbabilityDifferential, std::vector< highProbabilityDifferential >, std::function< bool(
        highProbabilityDifferential, highProbabilityDifferential ) > > pq( lambda );

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

/* Function that tries to find the best possible differential characteristic based on their probabilities of occurrence. */
void round_func_modified( uint32_t deltaX,
                          int roundNbr,
                          const std::vector< std::vector< std::vector< std::pair< uint32_t, int > >> >& sBoxDiffTables,
                          double currentProbability, std::vector< std::pair< uint32_t, uint32_t > > deltaXs )
{
    // Base case
    if ( roundNbr == FINAL_ROUND )
    {
        highProbabilityDifferential diff = { currentProbability, deltaXs };
        pq.push( diff );
        return;
    }

    // Extract the bits of all the s-boxes and determine which one is active
    uint32_t bitmask = 0xffffffff;

    uint32_t sBox1 = ( deltaX & ( bitmask << 28 ) ) >> 28;
    uint32_t sBox2 = ( deltaX & ( ( bitmask << 24 ) & ( bitmask >> 4 ) ) ) >> 24;
    uint32_t sBox3 = ( deltaX & ( ( bitmask << 20 ) & ( bitmask >> 8 ) ) ) >> 20;
    uint32_t sBox4 = ( deltaX & ( ( bitmask << 16 ) & ( bitmask >> 12 ) ) ) >> 16;
    uint32_t sBox5 = ( deltaX & ( ( bitmask << 12 ) & ( bitmask >> 16 ) ) ) >> 12;
    uint32_t sBox6 = ( deltaX & ( ( bitmask << 8 ) & ( bitmask >> 20 ) ) ) >> 8;
    uint32_t sBox7 = ( deltaX & ( ( bitmask << 4 ) & ( bitmask >> 24 ) ) ) >> 4;
    uint32_t sBox8 = deltaX & ( bitmask >> 28 );

    std::vector< std::pair< uint32_t, int > > sBox1DeltaY = sBoxDiffTables[ 0 ][ sBox1 ];
    std::vector< std::pair< uint32_t, int > > sBox2DeltaY = sBoxDiffTables[ 1 ][ sBox2 ];
    std::vector< std::pair< uint32_t, int > > sBox3DeltaY = sBoxDiffTables[ 2 ][ sBox3 ];
    std::vector< std::pair< uint32_t, int > > sBox4DeltaY = sBoxDiffTables[ 3 ][ sBox4 ];
    std::vector< std::pair< uint32_t, int > > sBox5DeltaY = sBoxDiffTables[ 4 ][ sBox5 ];
    std::vector< std::pair< uint32_t, int > > sBox6DeltaY = sBoxDiffTables[ 5 ][ sBox6 ];
    std::vector< std::pair< uint32_t, int > > sBox7DeltaY = sBoxDiffTables[ 6 ][ sBox7 ];
    std::vector< std::pair< uint32_t, int > > sBox8DeltaY = sBoxDiffTables[ 7 ][ sBox8 ];

    // Determine the number of active s-boxes
    std::vector< std::pair< uint8_t, std::vector< std::pair< uint32_t, int > > > > activeSBoxes;
    // the first parameter of the pair is the position of the s-box in the round function:
    // i.e. starting from left to right: s-box 1 is at position 7, s-box 2 is at position 6 ...
    // and s-box 8 is at position 0
    // This will help us determine how much we should shift the bits when calculating the result

    // Determine inactive s-boxes
    std::vector< uint8_t > inactiveSBoxes;

    if ( sBox1 > 0 )
        activeSBoxes.emplace_back( 7, sBox1DeltaY );
    else inactiveSBoxes.push_back( 7 );
    if ( sBox2 > 0 )
        activeSBoxes.emplace_back( 6, sBox2DeltaY );
    else inactiveSBoxes.push_back( 6 );
    if ( sBox3 > 0 )
        activeSBoxes.emplace_back( 5, sBox3DeltaY );
    else inactiveSBoxes.push_back( 5 );
    if ( sBox4 > 0 )
        activeSBoxes.emplace_back( 4, sBox4DeltaY );
    else inactiveSBoxes.push_back( 4 );
    if ( sBox5 > 0 )
        activeSBoxes.emplace_back( 3, sBox5DeltaY );
    else inactiveSBoxes.push_back( 3 );
    if ( sBox6 > 0 )
        activeSBoxes.emplace_back( 2, sBox6DeltaY );
    else inactiveSBoxes.push_back( 2 );
    if ( sBox7 > 0 )
        activeSBoxes.emplace_back( 1, sBox7DeltaY );
    else inactiveSBoxes.push_back( 1 );
    if ( sBox8 > 0 )
        activeSBoxes.emplace_back( 0, sBox8DeltaY );
    else inactiveSBoxes.push_back( 0 );

    // Generate all possibilities of active s-boxes to find the highest probability

    // vector of iterators to generate all possibilities
    std::vector< std::vector< std::pair< uint32_t, int > >::iterator > iterators( activeSBoxes.size() );

    // Instantiate all iterators
    for ( int i = 0; i < activeSBoxes.size(); ++i )
    {
        iterators[ i ] = activeSBoxes[ i ].second.begin();
    }

    // Odometer algorithm: https://stackoverflow.com/questions/1700079/howto-create-combinations-of-several-vectors-without-hardcoding-loops-in-c
    while ( iterators[ 0 ] != activeSBoxes[ 0 ].second.end() )
    {
        std::vector< double > probabilities;
        for ( const auto& it : iterators )
        {
            probabilities.push_back( it->second / 16.0 );
        }

        // Compute the probability of obtaining this result
        double updatedProbability = std::accumulate( probabilities.begin(), probabilities.end(), 1.0,
                                                     std::multiplies<>() ) * currentProbability;

        if ( updatedProbability > 0 )
        {
            // All probabilities are greater than 0. This means that we can generate a result
            uint32_t result = 0x00000000; // all non-active s-boxes will have zero at their position
            for ( int i = 0; i < iterators.size(); ++i )
            {
                uint32_t deltaY = iterators[ i ]->first;
                int shiftAmount = activeSBoxes[ i ].first * 4;
                result |= deltaY << shiftAmount;
            }

            // Now we can permute and recurse on this result
            uint32_t permutationResult = permute( result );

            deltaXs.emplace_back( deltaX, permutationResult );

            // The permutation result is going to be the input to the next round of the cipher
            round_func_modified( permutationResult, roundNbr + 1, sBoxDiffTables, updatedProbability, deltaXs );

            deltaXs.pop_back(); // backtrack
        }

        // Increment the odometer
        int K = ( int ) iterators.size();
        ++iterators[ K - 1 ]; // increment the last iterator
        for ( int i = K - 1; ( i > 0 ) && ( iterators[ i ] == activeSBoxes[ i ].second.end() ); --i )
        {
            iterators[ i ] = activeSBoxes[ i ].second.begin(); // reset the iterator
            ++iterators[ i - 1 ]; // increment the previous one
        }
    }
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

void crackKey4()
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

    // workerThread( 0, 0xffffffff, deltaPairs, lastRoundInputDifference );

    // Start first thread
    std::packaged_task< std::set< uint32_t >( uint32_t,
                                              uint32_t,
                                              std::vector< std::vector< uint64_t > >,
                                              uint32_t
    ) > task1 { workerThreadToCrackKey4 }; // Create task
    auto future1 = task1.get_future();      // Get the future object.

    // std::packaged_task is move-only
    std::thread t1( std::move( task1 ), 0, 0x40000000, deltaPairs, lastRoundInputDifference );

    // Start second thread
    std::packaged_task< std::set< uint32_t >( uint32_t,
                                              uint32_t,
                                              std::vector< std::vector< uint64_t > >,
                                              uint32_t
    ) > task2 { workerThreadToCrackKey4 }; // Create task
    auto future2 = task2.get_future();      // Get the future object.

    // std::packaged_task is move-only
    std::thread t2( std::move( task2 ), 0x40000000, 0x7FFFFFFF, deltaPairs, lastRoundInputDifference );

    // Start third thread
    std::packaged_task< std::set< uint32_t >( uint32_t,
                                              uint32_t,
                                              std::vector< std::vector< uint64_t > >,
                                              uint32_t
    ) > task3 { workerThreadToCrackKey4 }; // Create task
    auto future3 = task3.get_future();      // Get the future object.

    // std::packaged_task is move-only
    std::thread t3( std::move( task3 ), 0x7FFFFFFF, 0xBFFFFFFF, deltaPairs, lastRoundInputDifference );

    // Start fourth threads
    std::packaged_task< std::set< uint32_t >( uint32_t,
                                              uint32_t,
                                              std::vector< std::vector< uint64_t > >,
                                              uint32_t
    ) > task4 { workerThreadToCrackKey4 }; // Create task
    auto future4 = task4.get_future();      // Get the future object.

    // std::packaged_task is move-only
    std::thread t4( std::move( task4 ), 0xBFFFFFFF, 0xFFFFFFFF, deltaPairs, lastRoundInputDifference );

    t1.join();
    t2.join();
    t3.join();
    t4.join();

    std::set< uint32_t > s1 = future1.get();
    std::set< uint32_t > s2 = future2.get();
    std::set< uint32_t > s3 = future3.get();
    std::set< uint32_t > s4 = future4.get();

    std::set< uint32_t > subKeysK4;

    std::set_union( s1.begin(), s1.end(),
                    s2.begin(), s2.end(),
                    std::inserter( subKeysK4, std::begin( subKeysK4 ) )
    );

    std::set_union( s3.begin(), s3.end(),
                    s4.begin(), s4.end(),
                    std::inserter( subKeysK4, std::begin( subKeysK4 ) )
    );

    std::cout << "Final set has " << subKeysK4.size() << " sub-keys\n";

    // Write all sub-keys to a file
    std::ofstream file( "subKeysK4" );
    for ( uint32_t key : subKeysK4 )
    {
        file << std::hex << key << "\n";
    }
    file.close();
}

void workerThreadToCrackKey3(
        uint32_t begin,
        uint32_t end,
        std::set< uint32_t > subKeysK4,
        std::vector< std::vector< uint64_t > > pairs,
        uint32_t diffCharacteristicToRound3
)
{
    for ( uint32_t k3 = begin; k3 < end; ++k3 )
    {
        for ( uint32_t k4 : subKeysK4 )
        {
            int score = 0;
            for ( const std::vector< uint64_t >& pair : pairs )
            {
                // uint32_t key = 0x0BCA3B98;

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

                // uint32_t k3 = 0xA6C91854;
                uint32_t round3Result1 = round_func( result1, k3 );
                uint32_t round3Result2 = round_func( result2, k3 );

                uint32_t out = round3Result1 ^ round3Result2 ^ rightHalf;

                if ( diffCharacteristicToRound3 == out )
                {
                    // Before printing that we have found a match for k3, we need to make sure that bits in k3 match the
                    // corresponding bits of k4
                    bool allBitsMatch = true;
                    int bitPosition = 30;
                    for ( int i = 17; i <= 32; ++i )
                    {
                        // Extract the current bit of sub-key k3 and place it at the right position in the final integer
                        uint8_t ithBitOfK3 = ( k3 & ( 0x80000000 >> ( i - 1 ) ) ) > 0 ? 1 : 0;
                        uint8_t correspondingBitInK4 = ( k4 & ( 1 << bitPosition ) ) > 0 ? 1 : 0;

                        if ( ithBitOfK3 != correspondingBitInK4 )
                        {
                            allBitsMatch = false;
                            break;
                        }
                        bitPosition -= 2;
                    }

                    if ( allBitsMatch )
                    {
                        ++score;
                    }
                }
            }
            if ( score == pairs.size() )
            {
                std::cout << "Found sub-keys k3 & k4: " << "0x" << k3 << ", 0x" << k4 << std::endl;
            }
        }
    }
}

void crackKey3()
{
    // Let's break sub-key k3
    std::set< uint32_t > subKeysK4;
    std::ifstream file( "subKeysK4" );
    std::string subKey;
    while ( getline( file, subKey ) )
    {
        subKeysK4.insert( std::stoul( subKey, nullptr, 16 ) );
    }
    file.close();

    std::cout << "Finished loading sub-keys into set...proceeding with cracking round 3\n";

    // [Plaintext 1, Ciphertext 1, P2, C2]
    std::vector< std::vector< uint64_t > > pairs = {
            { 0xe0a621281e76a643, 0xFF3CCAC4EE3FFC8D, 0x1f59ded71e76a643, 0xB67A34C677C08000 },
            { 0x9d5da5a9eff0f367, 0xF735CA49CF34C3BD, 0x62a25a56eff0f367, 0xA2F7254BC9D2FEF9 },
            { 0x1629032bf2d0f365, 0x19CCF3B531745F36, 0xe9d6fcd4f2d0f365, 0x50BB3F6659250D1C },
            { 0xf2265d81fc1bd964, 0xBE9D1520BEA28942, 0x0dd9a27efc1bd964, 0x41FFD51AA46DBFB4 },
            { 0x03a54c0a20e8f05d, 0x59813499F0F5B72C, 0xfc5ab3f520e8f05d, 0xECF7D87793A6DB41 },
            { 0x237f61071881097b, 0xCA2FD612D8659769, 0xdc809ef81881097b, 0x86270DF629C2CD33 },
            { 0x44a3c3ecbdb612a2, 0xD6F435EF0712C355, 0xbb5c3c13bdb612a2, 0xD25C06AB2DC7F6AE },
            { 0x1f27d48d2c0fe9c8, 0x6EF64F5918F4138C, 0xe0d82b722c0fe9c8, 0x49C37771E5E57D8E },
            { 0x1196d39c49a80274, 0xBE4F6B195A65FE95, 0xee692c6349a80274, 0x3E9E2593047B520E },
            { 0x809606499db89b35, 0x2B90DDAF8BB4E8E7, 0x7f69f9b69db89b35, 0x0FF0939D84E7BF9A },
            { 0x81cbfdba3388b1ba, 0xF6053272C01CB9CB, 0x7e3402453388b1ba, 0x70AF7CAA283A124C },
            { 0x4a28d313c14be57e, 0x0A1EA4A827794303, 0xb5d72cecc14be57e, 0x84EE886538EA684D }
    };

    uint32_t diffCharacteristicToRound3 = 0xffffffff;
    std::cout << std::hex;

    std::thread t1( workerThreadToCrackKey3, 0, 0x40000000, subKeysK4, pairs, diffCharacteristicToRound3 );
    std::thread t2( workerThreadToCrackKey3, 0x40000000, 0x7FFFFFFF, subKeysK4, pairs, diffCharacteristicToRound3 );
    std::thread t3( workerThreadToCrackKey3, 0x7FFFFFFF, 0xBFFFFFFF, subKeysK4, pairs, diffCharacteristicToRound3 );
    std::thread t4( workerThreadToCrackKey3, 0xBFFFFFFF, 0xFFFFFFFF, subKeysK4, pairs, diffCharacteristicToRound3 );

    std::cout << "All threads have started...\n";

    t1.join();
    t2.join();
    t3.join();
    t4.join();
}

int main()
{
    // Part 1: Compute difference distribution tables for all s-boxes
    precompute_wes_permutation_mask(); /* Just an optimization: makes permutation step a bit faster */

    std::vector< std::vector< std::vector< int>> > diffTables = sBoxDifferenceDistributionTables();
    std::vector< std::vector< std::vector< std::pair< uint32_t, int > >> > optimizedDiffTables = optimizeDifferenceDistributionTables(
            diffTables );

    crackKey4();
    crackKey3();

    return 0;
}

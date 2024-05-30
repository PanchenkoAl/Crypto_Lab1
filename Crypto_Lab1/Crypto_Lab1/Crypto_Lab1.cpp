// Crypto_Lab1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <bitset>
#include <math.h>
#include <sstream>
#include <vector>
#include <cassert>

void FillTheMask(int* mask, int pow, int size);
long long int fastModulePower(int num, int pow, int module);
bool MillerRabin(int n);
int jacobi(int a, int b);
bool BPSW(int n);
void PrimeOutput(long long int n);
void runTests();

int gcd(int a, int b) 
{
    return a ? gcd(b % a, a) : b;
}

int main()
{
    long long p = 10;
    while(!BPSW(p))
        p = 10000 + (rand() % 10000);

    std::cout << fastModulePower(3, 5, 13) << std::endl;

    PrimeOutput(p);

    if (MillerRabin(23))
        std::cout << "is prime\n";
    else
        std::cout << "is not prime\n";
    if (BPSW(23))
        std::cout << "is prime\n";
    else
        std::cout << "is not prime\n";

    runTests();
}

void runTests()
{
    // Test gcd
    assert(gcd(48, 18) == 6);
    assert(gcd(101, 10) == 1);

    // Test fastModulePower
    assert(fastModulePower(5, 9, 1000) == 123);
    assert(fastModulePower(3, 5, 13) == 9);

    // Test jacobi
    assert(jacobi(1001, 9907) == -1);
    assert(jacobi(19, 45) == 1);

    // Test MillerRabin
    assert(MillerRabin(23) == true);
    assert(MillerRabin(24) == false);

    // Test BPSW
    assert(BPSW(23) == true);
    assert(BPSW(24) == false);

    // Test PrimeOutput function
    std::stringstream buffer;
    std::streambuf* old = std::cout.rdbuf(buffer.rdbuf());
    PrimeOutput(23);
    std::cout.rdbuf(old);

    std::string output = buffer.str();
    assert(output.find("base10: 23") != std::string::npos);
    assert(output.find("base2: 00000000000000000000000000010111") != std::string::npos);
    assert(output.find("base64: ") != std::string::npos);
    assert(output.find("byte[]: 23") != std::string::npos);
}

void PrimeOutput(long long int n)
{
    std::cout << "base10: " << n << std::endl;
    std::cout << "base2: " << std::bitset<32>(n) << std::endl;
    std::string base64;
    std::stringstream ss;
    ss << std::bitset<32>(n);
    std::bitset<6> bits;
    while (ss.good())
    {
        ss >> bits;
        int decimal = bits.to_ulong();
        base64 += "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[decimal];
    }
    std::cout << "base64: " << base64 << std::endl;
    std::vector<unsigned char> bytes;
    while (n > 0)
    {
        bytes.insert(bytes.begin(), n & 0xFF);
        n >>= 8;
    }
    std::cout << "byte[]: ";
    for (auto byte : bytes)
    {
        std::cout << static_cast<int>(byte) << " ";
    }
    std::cout << std::endl;
}

long long generatePrimeWithBits(int bits)
{
    long long lower = 1LL << (bits - 1);
    long long upper = (1LL << bits) - 1;

    long long candidate;
    do 
    {
        candidate = lower + (std::rand() % (upper - lower + 1));
        if (candidate % 2 == 0) candidate++;
    } while (!BPSW(candidate));

    return candidate;
}

bool BPSW(int n)
{
    if ((int)sqrt(n + 0.0) * (int)sqrt(n + 0.0) == n) return false;
    int dd = 5;
    for (;;)
    {
        int g = gcd(n, abs(dd));
        if (1 < g && g < n) return false;
        if (jacobi(dd, n) == -1) break;
        dd = dd < 0 ? -dd + 2 : -dd - 2;
    }
    int p = 1, q = (p * p - dd) / 4;
    int d = n + 1, s = 0;
    while ((d & 1) == 0)
        ++s, d >>= 1;
    long long u = 1, v = p, u2m = 1, v2m = p, qm = q, qm2 = q * 2, qkd = q;
    for (int mask = 2; mask <= d; mask <<= 1)
    {
        u2m = (u2m * v2m) % n;
        v2m = (v2m * v2m) % n;
        while (v2m < qm2)   v2m += n;
        v2m -= qm2;
        qm = (qm * qm) % n;
        qm2 = qm * 2;
        if (d & mask) 
        {
            long long t1 = (u2m * v) % n, t2 = (v2m * u) % n,
                t3 = (v2m * v) % n, t4 = (((u2m * u) % n) * dd) % n;
            u = t1 + t2;
            if (u & 1)  u += n;
            u = (u >> 1) % n;
            v = t3 + t4;
            if (v & 1)  v += n;
            v = (v >> 1) % n;
            qkd = (qkd * qm) % n;
        }
    }
    if (u == 0 || v == 0)  return true;
    long long qkd2 = qkd * 2;
    for (int r = 1; r < s; ++r) {
        v = (v * v) % n - qkd2;
        if (v < 0)  v += n;
        if (v < 0)  v += n;
        if (v >= n)  v -= n;
        if (v >= n)  v -= n;
        if (v == 0)  return true;
        if (r < s - 1) {
            qkd = (qkd * 1ll * qkd) % n;
            qkd2 = qkd * 2;
        }
    }
    return false;
}

int jacobi(int a, int b)
{
    if (a == 0) return 0;
    if (a == 1) return 1;
    if (a < 0)
        if ((b & 2) == 0)
            return jacobi(-a, b);
        else
            return -jacobi(-a, b);
    int a1 = a, e = 0;
    while ((a1 & 1) == 0)
        a1 >>= 1, ++e;
    int s;
    if ((e & 1) == 0 || (b & 7) == 1 || (b & 7) == 7)
        s = 1;
    else
        s = -1;
    if ((b & 3) == 3 && (a1 & 3) == 3)
        s = -s;
    if (a1 == 1)
        return s;
    return s * jacobi(b % a1, a1);
}

bool MillerRabin(int n)
{
    if ((n <= 2) || (n % 2 == 0))
    {
        std::cout << "Wrong input, try again.\n";
        return 0;
    }

    int s = 1;
    int t = n;
    for (;;s++)
    {
        t /= 2;
        if (t % 2 == 1)
            break;
    }

    for (int i = 0; i < log2(n)+1; i++)
    {
        int a = rand() % n;

        int x = fastModulePower(a, t, n);
        if ((x == 1) || (x == n - 1))
            continue;
        for (int j = 0; j < s - 1; j++)
        {
            x = fastModulePower(x, 2, n);
            if (x == 1)
                return false;
            if (x == n - 1)
                break;
        }
        return false;
    }
    return true;
}

long long int fastModulePower(int num, int pow, int module)
{
    int mask_size = log2(pow) + 1;
    int* mask = new int[mask_size];
    int* modules = new int[mask_size];
    long long int result = 1;

    FillTheMask(mask, pow, mask_size);

    modules[0] = num;
    for (int i = 1; i < mask_size; i++)
    {
        modules[i] = (modules[i - 1] * modules[i - 1]) % module;
    }

    for (int i = 0; i < mask_size; i++)
    {
        if (mask[i] == 1)
            result *= modules[i];
    }

    return result % module;
}

void FillTheMask(int* mask, int pow, int size)
{
    int check;
    double d = pow;
    for (int i = 0; i < size; i++)
    {
        check = 1;
        for (int j = 0; j < size; j++)
        {
            if (check > d/2)
            {
                mask[j] = 1;
                d -= check;
                break;
            }
            check *= 2;
        }
    }
}
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

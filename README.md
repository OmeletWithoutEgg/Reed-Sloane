# Reed-Sloane Algorithm
keywords: LFSR, Berlekamp-Massey, linear recurrence

[Reed-Sloane Algorithm](https://en.wikipedia.org/wiki/Reeds%E2%80%93Sloane_algorithm) is an extension to the Berlekamp-Massey algorithm which applies when the terms of the sequences are integers modulo some given modulus m.
The implementation is almost based on [this paper](http://neilsloane.com/doc/Me111.pdf).
The time complexity is O(Ω(m)n^2), where Ω(m) is the total number of prime factors of m and n is the length of the given sequence.
There are some details in this implementation:
1. We don't have to calculate b(x). Instead, we maintain L(a(x), b(x)), since the exact value of b(x) doesn't matter when we calculating the x^k term of (S * a(x) - b(x)).
2. When calculating (S * a(x) - b(x)), we can take less modulo operations, see the code for more details.

You can solve [this problem](https://www.spoj.com/problems/FINDLR/)
by using Reed-Sloane Algorithm, but must be careful that modulo operations may overflow in such constraints.

# Open track: query-cost model compared with measurement

## Cost of one read

Reading a phase to L bits with semiclassical phase estimation uses powers 2^(L-1), …, 2, 1. Each gate costs
bitlen(power), so one read costs

    sum_{k=0}^{L-1} (k + 1) = L(L+1)/2.

Each oracle gate is one control qubit, so it yields at most one bit of information. A high-precision bit needs
a large power, so a read cannot be much cheaper than this quadratic cost.

## Happy-path prediction (b = bitlen(n))

| Step | Bits / shots | Queries |
| --- | --- | --- |
| Order read (g) | L = 126 | 8001 |
| Order test g^n = e | 16 shots | 16·b |
| Log read (g, A, B, C) | 4 reads of b+7 bits | 4·(b+7)(b+8)/2 |
| Three log tests | 16 shots each | 3·16·(b+1) |

For public_1 (b = 56) this gives 8001 + 896 + 8064 + 2736 = **19 697**. The checker measured **19 601**, the
most common value in `results.txt`. The difference is that the test exponents are reduced to (−n/2, n/2], so
they have fewer bits.

## Why measured costs are sometimes higher

- **Noise.** The phase error of a gate has standard deviation 0.002·√cost turns, and it corrupts that qubit
  with probability about π²σ². Summed over a 126-bit read, P(bad order read) ≈ 3.95e-5 · Σ(k+1) ≈ 0.3.
  A bad read only fails the identity test, and the solver tries again.
- **gcd(j, n) > 1.** This happens with probability 1 − φ(n)/n. For public_2, n = 3³·5³·7²·11³·13², so it
  happens 62% of the time and triggers the lift step or a second sample.
- **Forgery.** Three failed samples are required before calling a forgery. This keeps the chance of
  calling a valid element forged below about 0.1%.

## How the solver degrades with noise

Expected failures scale as noise² · L². Doubling the noise makes about 4× as many reads go bad. The answer
stays correct, because every result is checked, but the solver spends more of its budget on retries.

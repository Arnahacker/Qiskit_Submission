"""The Silent Group solver (QFF hackathon, Track 2).

Each shot the device samples one character chi of the hidden group G, and
oracle_gate(h, k) kicks the phase k * chi(h) (in turns) back onto its control.
Every oracle gate in one shot sees the same chi, so one shot can read chi(g),
chi(A), chi(B) and chi(C) together.

1. Order. A semiclassical phase estimation reads chi(g) = j/n to 126 bits.
   It uses one qubit, reset each round, with if_test phase corrections.
   Continued fractions then give d = n / gcd(j, n). If g^d != e, the missing
   factor gcd(j, n) is the order of g^d. That order is small, so a few cheap
   phase reads of g^d recover it.
2. Logs. Knowing n, read chi(g) = j/n and chi(X) = x*j/n (X = A, B, C) to
   bitlen(n) + 7 bits. Then x = k_X / k_g mod n. If s = gcd(k_g, n) > 1, x is
   only known mod n/s: brute-force a small s, or CRT with another sample.
3. Every answer is checked with an identity test: H, oracle(s), H, measure.
   If the element is e, the qubit always reads 0. Otherwise each shot reads 1
   with probability exactly 1/2. So 16 shots with at most one 1 leave a false
   accept chance of 17/2^16.
   An element is called a forgery only after three independent samples fail.
   This keeps noisy reads from looking like forgeries.
"""

from __future__ import annotations

import math

from qiskit import QuantumCircuit

from silentgroup import oracle_gate

MAX_ORDER = 1 << 60  # |G| <= 2^60, so ord(g) <= 2^60
ORDER_BITS = 2 * 60 + 6  # continued fractions need |error| < 1/(2 n^2); 6 spare bits
LOG_MARGIN = 7  # extra bits over bitlen(n) when reading chi(X) = k/n
CORRECTION_WINDOW = 12  # semiclassical QFT: drop corrections below 2^-13 turns
TEST_SHOTS = 16
TEST_MAX_ONES = 1  # noise flips an identity test shot with probability ~0.2%
LIFT_MAX = 1 << 12  # largest gcd(j, n) the cheap lift looks for
BRUTE_FORCE = 4  # try up to this many log candidates before taking another sample
FORGERY_FAILURES = 3
MAX_LOG_SAMPLES = 8


class _OutOfBudget(Exception):
    pass


def _run(oracle, qc, shots):
    left = oracle.budget_left()
    if left is not None and oracle.cost(qc, shots) > left:
        raise _OutOfBudget
    return oracle.run(qc, shots=shots)


# ---------------------------------------------------------------- circuits


def _append_oracles(qc, terms, scale):
    """Multiply by prod(elem ** (mult * scale)) controlled on qubit 0."""
    for elem, mult in terms:
        if mult * scale:
            qc.append(oracle_gate(elem, mult * scale), [0])


def _append_read(qc, terms, nbits, clbit0):
    """Semiclassical phase estimation of chi(prod elem^mult) into clbits clbit0.. (LSB first).

    Round i uses power 2^(nbits-1-i), whose phase is m/2^(i+1) mod 1. Measured bits
    0..i-1 cancel the fractional tail, so the qubit is left with phase 0 or 1/2.
    """
    for i in range(nbits):
        qc.h(0)
        _append_oracles(qc, terms, 1 << (nbits - 1 - i))
        for b in range(max(0, i - CORRECTION_WINDOW), i):
            with qc.if_test((qc.clbits[clbit0 + b], 1)):
                qc.p(-2 * math.pi / 2 ** (i + 1 - b), 0)
        qc.h(0)
        qc.measure(0, clbit0 + i)
        qc.reset(0)


def _read(oracle, reads, shots=1):
    """Run phase reads [(terms, nbits), ...] in one circuit (one chi per shot).

    Returns one list of integers m (phase ~ m / 2^nbits) per shot.
    """
    qc = QuantumCircuit(1, sum(nbits for _, nbits in reads))
    spans, offset = [], 0
    for terms, nbits in reads:
        _append_read(qc, terms, nbits, offset)
        spans.append((offset, nbits))
        offset += nbits
    out = []
    for key, count in _run(oracle, qc, shots).items():
        n = len(key)  # clbit 0 is the rightmost character
        out += [[int(key[n - o - nb : n - o], 2) for o, nb in spans]] * count
    return out


def _is_identity(oracle, terms):
    """Test prod(elem^mult) == e: P(read 1) is ~0 if so, and exactly 1/2 per shot if not."""
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    _append_oracles(qc, terms, 1)
    qc.h(0)
    qc.measure(0, 0)
    return _run(oracle, qc, TEST_SHOTS).get("1", 0) <= TEST_MAX_ONES


# ---------------------------------------------------------------- arithmetic


def _denominator(m, nbits, qmax):
    """Largest continued-fraction convergent denominator <= qmax of m / 2^nbits."""
    num, den = m, 1 << nbits
    q0, q1, best = 1, 0, 1
    while den:
        a = num // den
        num, den = den, num - a * den
        q0, q1 = q1, a * q1 + q0
        if q1 > qmax:
            break
        best = q1
    return best


def _lattice(m, nbits, n):
    """k in Z_n with m / 2^nbits ~ k / n, or None if the read is far from every k / n."""
    num = m * n
    k, rem = divmod(num + (1 << (nbits - 1)), 1 << nbits)
    off = rem - (1 << (nbits - 1))  # offset from k, in units of 2^-nbits / n
    return k % n if 4 * abs(off) <= (1 << nbits) else None


def _small_factors(s):
    """Prime factors of s by trial division; a large cofactor is returned as is."""
    out, p = [], 2
    while p * p <= s and p < 1 << 16:
        if s % p == 0:
            out.append(p)
            while s % p == 0:
                s //= p
        p += 1 if p == 2 else 2
    if s > 1:
        out.append(s)
    return out


def _crt(r1, m1, r2, m2):
    g = math.gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    lcm = m1 // g * m2
    t = (r2 - r1) // g * pow(m1 // g, -1, m2 // g) % (m2 // g)
    return (r1 + m1 * t) % lcm, lcm


def _sym(x, n):
    """Representative of x mod n with the smallest |x|, so oracle powers stay cheap."""
    x %= n
    return x - n if x > n // 2 else x


# ---------------------------------------------------------------- order


def _lift(oracle, d):
    """g^d has order t = gcd(j, n); read it from two cheap samples.

    t is almost always small, so read only enough bits for t <= LIFT_MAX; a larger t
    falls through to another full order read.
    """
    tmax = min(MAX_ORDER // d, LIFT_MAX)
    if tmax < 2:
        return None
    nbits = 2 * tmax.bit_length() + 5
    t = 1
    for (m,) in _read(oracle, [([("g", d)], nbits)], shots=2):
        t = math.lcm(t, _denominator(m, nbits, tmax))
    if 1 < t and d * t <= MAX_ORDER and _is_identity(oracle, [("g", d * t)]):
        return d * t
    return None


def find_order(oracle):
    tried, seen = set(), []
    while True:
        ((m,),) = _read(oracle, [([("g", 1)], ORDER_BITS)])
        d = _denominator(m, ORDER_BITS, MAX_ORDER)
        for c in [d] + [math.lcm(d, e) for e in seen[-2:]]:
            if c <= MAX_ORDER and c not in tried:
                tried.add(c)
                if _is_identity(oracle, [("g", c)]):
                    return c
        c = _lift(oracle, d)
        if c:
            return c
        seen.append(d)


# ---------------------------------------------------------------- logs


def _candidates(rows, n):
    """Residue classes for x from single samples and pairs, small enough to brute force."""
    out = []
    for i, (r, m) in enumerate(rows):
        out.append((r, m, i))
        for r2, m2 in rows[:i]:
            merged = _crt(r, m, r2, m2)
            if merged:
                out.append((*merged, None))
    return [(r, m, i) for r, m, i in out if n // m <= BRUTE_FORCE]


def find_logs(oracle, n):
    """Return (n, {X: log}) with X = g^log, or (n, None) if some X is a forgery.

    n may shrink: a sample can reveal that a lifted order still had a spare factor.
    """
    if n == 1:
        return n, ({X: 0 for X in "ABC"} if all(_is_identity(oracle, [(X, 1)]) for X in "ABC") else None)
    pending = ["A", "B", "C"]
    logs, rows = {}, {X: [] for X in pending}
    tested = {X: set() for X in pending}
    failures = dict.fromkeys(pending, 0)
    checked_primes = set()
    for _ in range(MAX_LOG_SAMPLES):
        nbits = n.bit_length() + LOG_MARGIN
        (vals,) = _read(oracle, [([("g", 1)], nbits)] + [([(X, 1)], nbits) for X in pending])
        kg = _lattice(vals[0], nbits, n)
        if kg is None:
            continue  # chi(g) misread: says nothing about A, B, C
        # If n were k * ord(g), every k_g would be a multiple of k: check primes of gcd(k_g, n).
        # A noisy order read can land on lcm(ord(g), 2^i), so divide out whole prime powers.
        shrunk = False
        for p in _small_factors(math.gcd(kg, n)):
            while n % p == 0 and p not in checked_primes:
                if _is_identity(oracle, [("g", n // p)]):
                    n //= p
                    shrunk = True
                else:
                    checked_primes.add(p)
        if shrunk:  # this sample was read with more than enough bits for the smaller n
            rows, tested = {X: [] for X in "ABC"}, {X: set() for X in "ABC"}
            kg = _lattice(vals[0], nbits, n)
            if kg is None:
                continue
        s = math.gcd(kg, n)
        inv = pow(kg // s, -1, n // s)
        for idx, X in enumerate(list(pending)):
            kx = _lattice(vals[idx + 1], nbits, n)
            if kx is None or kx % s:
                failures[X] += 1
            else:
                rows[X].append((kx // s * inv % (n // s), n // s))
                sample = len(rows[X]) - 1
                for r, m, origin in _candidates(rows[X], n):
                    xs = [r + m * u for u in range(n // m)]
                    found = None
                    for x in xs:
                        if x not in tested[X]:
                            tested[X].add(x)
                            if _is_identity(oracle, [(X, 1), ("g", _sym(-x, n))]):
                                found = x
                                break
                    if found is not None:
                        logs[X] = found
                        pending.remove(X)
                        break
                    if origin == sample:
                        failures[X] += 1
            if X in pending and failures[X] >= FORGERY_FAILURES:
                return n, None
        if not pending:
            return n, logs
    return n, None


def solve(oracle) -> dict:
    order = 1
    try:
        order = find_order(oracle)
        order, logs = find_logs(oracle, order)
    except _OutOfBudget:
        return {"order": order, "secret": "none"}
    if logs is None:
        return {"order": order, "secret": "none"}
    return {"order": order, "secret": logs["A"] * logs["B"] * logs["C"] % order}

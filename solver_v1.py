"""The Silent Group Solver.

Track 2 Submission for Qiskit Fall Fest 2026.
Uses multi-stage overlapping Quantum Phase Estimation (QPE) with classical backward stitching,
Shor continued fractions order extraction with smooth-multiplier search and multi-shot LCM,
modular discrete logarithm recovery, and black-box forgery verification.
"""

from __future__ import annotations

import math
from fractions import Fraction

from qiskit import QuantumCircuit
from qiskit.circuit.library import QFTGate

from silentgroup import oracle_gate

T_QUBITS = 14
STEP_BITS = 6
G_STAGES = 19
LOG_STAGES = 11

SMOOTH_MULTIPLIERS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 24, 25, 27, 30, 32, 35, 36, 40, 42, 45, 48, 49, 50]


def quick_verify_order(oracle, n: int) -> bool:
    """Check whether g^n == e on the group register with 2-stage verification."""
    if n <= 1000 or n > 2**60:
        return False
    b = oracle.budget_left()
    if b is not None and b < 900:
        return False

    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.append(oracle_gate("g", n), [0])
    qc.h(0)
    qc.measure(0, 0)
    counts = oracle.run(qc, shots=3)
    if counts.get("0", 0) < 3:
        return False
    counts2 = oracle.run(qc, shots=10)
    return counts2.get("0", 0) >= 9


def verify_log(oracle, elem: str, exp: int) -> bool:
    """Verify whether elem == g^exp by testing (elem * g^(-exp)) == e with 2-stage filter."""
    b = oracle.budget_left()
    if b is not None and b < 800:
        return False

    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.append(oracle_gate(elem, 1), [0])
    qc.append(oracle_gate("g", -exp), [0])
    qc.h(0)
    qc.measure(0, 0)
    counts = oracle.run(qc, shots=2)
    if counts.get("0", 0) < 2:
        return False
    counts2 = oracle.run(qc, shots=9)
    return counts2.get("0", 0) >= 8


def round_fraction(frac: Fraction, order: int) -> int:
    val = (frac % 1) * order
    return (val.numerator + val.denominator // 2) // val.denominator


def reconstruct_phase(raw: str, start_c: int, ns: int) -> Fraction:
    stages = []
    for s in range(ns):
        c_start = start_c + s * T_QUBITS
        chunk = raw[len(raw) - (c_start + T_QUBITS) : len(raw) - c_start]
        stages.append(int(chunk, 2))

    phi = Fraction(stages[-1], 2**T_QUBITS)
    for s in reversed(range(ns - 1)):
        coarse = Fraction(stages[s], 2**T_QUBITS)
        best_k = 0
        best_diff = 1.0
        for k in range(2**STEP_BITS):
            cand = (k + phi) / (2**STEP_BITS)
            diff = float(abs((cand - coarse + Fraction(1, 2)) % 1 - Fraction(1, 2)))
            if diff < best_diff:
                best_diff = diff
                best_k = k
        phi = (best_k + phi) / (2**STEP_BITS)
    return phi


def get_g_candidate_denominators(oracle) -> list[int]:
    """Run QPE on g with 19 stages (~122 bits) and return top candidate denominators."""
    total_clbits = G_STAGES * T_QUBITS
    qc = QuantumCircuit(T_QUBITS, total_clbits)

    for s in range(G_STAGES):
        shift = s * STEP_BITS
        for j in range(T_QUBITS):
            qc.h(j)
            qc.append(oracle_gate("g", 2 ** (shift + j)), [j])
        qc.append(QFTGate(T_QUBITS).inverse(), range(T_QUBITS))
        for j in range(T_QUBITS):
            qc.measure(j, s * T_QUBITS + j)
        if s < G_STAGES - 1:
            qc.reset(range(T_QUBITS))

    counts = oracle.run(qc, shots=1)
    raw = list(counts.keys())[0]
    phi_g = reconstruct_phase(raw, 0, G_STAGES)

    p0, q0 = phi_g.numerator, phi_g.denominator
    h, k = [0, 1], [1, 0]
    entries = []
    while q0:
        a = p0 // q0
        p0, q0 = q0, p0 % q0
        h.append(a * h[-1] + h[-2])
        k.append(a * k[-1] + k[-2])
        entries.append((a, h[-1], k[-1]))

    cands = []
    for idx, (a, num, den) in enumerate(entries):
        if 1000 < den <= 2**60:
            next_a = entries[idx + 1][0] if idx + 1 < len(entries) else 1
            cands.append((next_a, den))
    cands.sort(key=lambda x: x[0], reverse=True)

    return [den for _, den in cands[:6]]


def find_order(oracle) -> int | None:
    dens1 = get_g_candidate_denominators(oracle)

    # Fast pass: m = 1 across top candidates
    for d in dens1:
        if quick_verify_order(oracle, d):
            return d

    # Multiplier pass: smooth multipliers on top 3 candidates
    for d in dens1[:3]:
        for m in SMOOTH_MULTIPLIERS[1:]:
            cand = d * m
            if 1 < cand <= 2**60:
                if quick_verify_order(oracle, cand):
                    return cand

    # Shot 2 if needed and budget permits
    budget = oracle.budget_left()
    if budget is not None and budget >= 35000:
        dens2 = get_g_candidate_denominators(oracle)
        for d in dens2:
            if quick_verify_order(oracle, d):
                return d

        for d2 in dens2[:3]:
            for m in SMOOTH_MULTIPLIERS[1:]:
                cand = d2 * m
                if 1 < cand <= 2**60:
                    if quick_verify_order(oracle, cand):
                        return cand

        for d1 in dens1[:3]:
            for d2 in dens2[:3]:
                L = math.lcm(d1, d2)
                if 1000 < L <= 2**60:
                    if quick_verify_order(oracle, L):
                        return L

    return None


def solve_discrete_logs(oracle, order: int) -> dict[str, int | None]:
    """Measure g, A, B, C in a single shot with 11 stages to extract logs."""
    elements = ["g", "A", "B", "C"]
    total_clbits = len(elements) * LOG_STAGES * T_QUBITS
    qc = QuantumCircuit(T_QUBITS, total_clbits)
    clbit_idx = 0
    ranges = {}

    for elem in elements:
        start_c = clbit_idx
        for s in range(LOG_STAGES):
            shift = s * STEP_BITS
            for j in range(T_QUBITS):
                qc.h(j)
                qc.append(oracle_gate(elem, 2 ** (shift + j)), [j])
            qc.append(QFTGate(T_QUBITS).inverse(), range(T_QUBITS))
            for j in range(T_QUBITS):
                qc.measure(j, clbit_idx + j)
            clbit_idx += T_QUBITS
            qc.reset(range(T_QUBITS))
        ranges[elem] = start_c

    counts = oracle.run(qc, shots=1)
    raw = list(counts.keys())[0]

    phases = {elem: reconstruct_phase(raw, ranges[elem], LOG_STAGES) for elem in elements}

    base_kg = round_fraction(phases["g"], order)
    base_kes = {elem: round_fraction(phases[elem], order) for elem in ["A", "B", "C"]}

    best_inv_kg = None
    known_logs = {}

    # Anchor search: check delta_g = 0 first across all elements, then +-1, +-2
    for delta_g in (0, 1, -1, 2, -2):
        kg = (base_kg + delta_g) % order
        if math.gcd(kg, order) != 1:
            continue
        inv_kg = pow(kg, -1, order)

        for elem in ["A", "B", "C"]:
            ke = base_kes[elem]
            for delta_e in (0, 1, -1, 2, -2, 3, -3):
                cand = ((ke + delta_e) * inv_kg) % order
                if verify_log(oracle, elem, cand):
                    best_inv_kg = inv_kg
                    known_logs[elem] = cand
                    break
            if best_inv_kg is not None:
                break
        if best_inv_kg is not None:
            break

    if best_inv_kg is None:
        if math.gcd(base_kg, order) == 1:
            best_inv_kg = pow(base_kg, -1, order)
        else:
            return {"A": None, "B": None, "C": None}

    # Resolve remaining elements
    logs = {}
    for elem in ["A", "B", "C"]:
        if elem in known_logs:
            logs[elem] = known_logs[elem]
            continue
        ke = base_kes[elem]
        cand_log = None
        for delta_e in (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6):
            cand = ((ke + delta_e) * best_inv_kg) % order
            if verify_log(oracle, elem, cand):
                cand_log = cand
                break
        logs[elem] = cand_log

    return logs


def solve(oracle) -> dict:
    order = find_order(oracle)
    if order is None:
        return {"order": 1, "secret": "none"}

    b = oracle.budget_left()
    if b is not None and b < 24000:
        return {"order": order, "secret": "none"}

    logs = solve_discrete_logs(oracle, order)

    if any(logs[e] is None for e in ["A", "B", "C"]):
        secret = "none"
    else:
        secret = (logs["A"] * logs["B"] * logs["C"]) % order

    return {"order": order, "secret": secret}

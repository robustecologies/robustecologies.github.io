/-
  Exact arithmetic of the erratum to Almaraz, Kalita, Langa and Soler-Toscano (2024),
  Journal of Mathematical Biology 88, 64, checked by the kernel of Lean 4.35.0.

  Scope. Every theorem below is closed by `decide`, so the Lean kernel itself evaluates
  the claim; no compiled code and no external library is trusted. The file certifies
  the exact arithmetic on which the corrections rest: equilibria, the signs of the
  invasion rates, the admissible communities and the edges of the invasion graphs, the
  matrix counterexample to Lemma 11, the Sylvester minors of the Volterra-Lyapunov
  certificates, the residual of the vector printed with Figure 3, the rates, certificate,
  invasion graph and sign graph of connections of the threshold example of Question 27, and
  the invasion graph, the Hurwitz resident blocks and the invasion rate r2({1,3,4}) of the
  D-stable example of Question 35. The dynamical arguments of the erratum (attraction,
  connecting orbits, the proofs for Theorems 18, 23, 28 and 32, Lemma 37) are proved in the
  text and are not formalised here; that the sign graph of Question 27 is the graph of
  connections is Theorem 3.2 of Bortolan, Kalita, Langa and Moura (2026). The D-stability of the Question 35 matrix, the existence of its
  periodic orbit and the missing connection rest on three certificates in exact rational
  and ball arithmetic, which an independent verifier kept with the companion project
  checks outside this file and outside the reproducibility archive.

  Conventions. Species are indexed from 0, so the community {1,3} of the erratum is
  [0, 2] here. A matrix A and a growth vector b are multiplied by one positive integer
  so that every entry is an integer; this leaves every equilibrium unchanged and
  multiplies every invasion rate by that integer, so supports, signs and edges are
  those of the printed system. For a support I the equilibrium is obtained by Cramer's
  rule, x_k = N_k / D with D = det A(I), and an invasion rate b_i + sum_k a_{i,I_k} x_k is
  the quotient (b_i D + sum_k a_{i,I_k} N_k) / D, whose sign is that of the product of
  numerator and D.

  Run:  lean ErrataChecks.lean   (no output means that every theorem was accepted)
-/

abbrev Mat := List (List Int)
abbrev Vec := List Int

def vget (v : Vec) (i : Nat) : Int := v.getD i 0
def entry (A : Mat) (i j : Nat) : Int := vget (A.getD i []) j

def dropAt {T : Type} : List T -> Nat -> List T
  | [], _ => []
  | _ :: t, 0 => t
  | h :: t, k + 1 => h :: dropAt t k

def setAt : Vec -> Nat -> Int -> Vec
  | [], _, _ => []
  | _ :: t, 0, x => x :: t
  | h :: t, k + 1, x => h :: setAt t k x

/-- Determinant by Laplace expansion along the first row. -/
def detAux : Nat -> Mat -> Int
  | 0, _ => 1
  | n + 1, A =>
    match A with
    | [] => 1
    | r :: rest =>
      (List.range (n + 1)).foldl
        (fun acc j => acc + (if j % 2 == 0 then 1 else -1) * vget r j * detAux n (rest.map (fun s => dropAt s j))) 0

def det (A : Mat) : Int := detAux A.length A

def sub (A : Mat) (I : List Nat) : Mat := I.map (fun i => I.map (fun j => entry A i j))

def replCol (M : Mat) (k : Nat) (v : Vec) : Mat := (M.zip v).map (fun p => setAt p.1 k p.2)

def cramerD (A : Mat) (I : List Nat) : Int := det (sub A I)

def cramerN (A : Mat) (b : Vec) (I : List Nat) (k : Nat) : Int :=
  det (replCol (sub A I) k (I.map (fun i => - vget b i)))

/-- The support I carries a positive equilibrium. -/
def feasible (A : Mat) (b : Vec) (I : List Nat) : Bool :=
  let D := cramerD A I
  D != 0 && (List.range I.length).all (fun k => cramerN A b I k * D > 0)

/-- Numerator of the invasion rate of species i at the equilibrium of I (denominator D). -/
def rateNum (A : Mat) (b : Vec) (I : List Nat) (i : Nat) : Int :=
  vget b i * cramerD A I + (List.range I.length).foldl (fun acc k => acc + entry A i (I.getD k 0) * cramerN A b I k) 0

/-- A number with the sign of the invasion rate of species i at I. -/
def rateSign (A : Mat) (b : Vec) (I : List Nat) (i : Nat) : Int := rateNum A b I i * cramerD A I

def combos : Nat -> List Nat -> List (List Nat)
  | 0, _ => [[]]
  | _ + 1, [] => []
  | k + 1, h :: t => (combos k t).map (fun l => h :: l) ++ combos (k + 1) t

def allSubsets (n : Nat) : List (List Nat) := (List.range (n + 1)).flatMap (fun k => combos k (List.range n))

def admissible (A : Mat) (b : Vec) (n : Nat) : List (List Nat) := (allSubsets n).filter (feasible A b)

def absent (n : Nat) (I : List Nat) : List Nat := (List.range n).filter (fun i => !I.contains i)

def isEdge (A : Mat) (b : Vec) (I K : List Nat) : Bool :=
  I != K &&
  (K.filter (fun j => !I.contains j)).all (fun j => rateSign A b I j > 0) &&
  (I.filter (fun j => !K.contains j)).all (fun j => rateSign A b K j < 0)

def edges (A : Mat) (b : Vec) (n : Nat) : List (Prod (List Nat) (List Nat)) :=
  let E := admissible A b n
  E.flatMap (fun I => (E.filter (fun K => isEdge A b I K)).map (fun K => (I, K)))

/-- Sign condition of the graph of connections under Volterra-Lyapunov stability (Theorem 3.2
    of Bortolan, Kalita, Langa and Moura, 2026): positive entry rates at the source and
    nonpositive exclusion rates at the target. -/
def isConnEdge (A : Mat) (b : Vec) (I K : List Nat) : Bool :=
  I != K &&
  (K.filter (fun j => !I.contains j)).all (fun j => rateSign A b I j > 0) &&
  (I.filter (fun j => !K.contains j)).all (fun j => rateSign A b K j <= 0)

def connEdges (A : Mat) (b : Vec) (n : Nat) : List (Prod (List Nat) (List Nat)) :=
  let E := admissible A b n
  E.flatMap (fun I => (E.filter (fun K => isConnEdge A b I K)).map (fun K => (I, K)))

/-- Kahn's test: remove the vertices without incoming edges until none is left (acyclic) or
    every remaining vertex has an incoming edge (a cycle exists). Each round removes at least
    one vertex, so the number of vertices bounds the number of rounds. -/
def acyclicAux : Nat -> List (List Nat) -> List (Prod (List Nat) (List Nat)) -> Bool
  | 0, V, _ => V.isEmpty
  | k + 1, V, E =>
    let S := V.filter (fun v => !(E.any (fun e => e.2 == v)))
    if S.isEmpty then V.isEmpty
    else acyclicAux k (V.filter (fun v => !S.contains v)) (E.filter (fun e => !S.contains e.1))

def acyclic (V : List (List Nat)) (E : List (Prod (List Nat) (List Nat))) : Bool := acyclicAux V.length V E

/-- Every missing-species invasion rate at every admissible community is nonzero. -/
def nonzeroRates (A : Mat) (b : Vec) (n : Nat) : Bool :=
  (admissible A b n).all (fun I => (absent n I).all (fun i => rateSign A b I i != 0))

/-- Invasion index: the number of absent species with a positive invasion rate. -/
def index (A : Mat) (b : Vec) (n : Nat) (I : List Nat) : Nat :=
  ((absent n I).filter (fun i => rateSign A b I i > 0)).length

def negSym (A : Mat) (n : Nat) : Mat :=
  (List.range n).map (fun i => (List.range n).map (fun j => - (entry A i j + entry A j i)))

/-- Leading principal minors of a square matrix of size n. -/
def leadingMinors (M : Mat) (n : Nat) : List Int := (List.range n).map (fun k => det (sub M (List.range (k + 1))))

def shift (M : Mat) (n : Nat) (c : Int) : Mat :=
  (List.range n).map (fun i => (List.range n).map (fun j => entry M i j - (if i == j then c else 0)))

/-- Sum of the principal minors of order j of a square matrix. -/
def minorSum (M : Mat) (j : Nat) : Int :=
  (combos j (List.range M.length)).foldl (fun acc I => acc + det (sub M I)) 0

/-- Coefficient c_j of the characteristic polynomial l^k + c_1 l^(k-1) + ... + c_k. -/
def charCoeff (M : Mat) (j : Nat) : Int := (if j % 2 == 0 then 1 else -1) * minorSum M j

/-- Lienard-Chipart test: every eigenvalue of M (of size 1 to 4) has negative real part. -/
def hurwitz (M : Mat) : Bool :=
  let c := charCoeff M
  match M.length with
  | 1 => c 1 > 0
  | 2 => c 1 > 0 && c 2 > 0
  | 3 => c 1 > 0 && c 3 > 0 && c 1 * c 2 > c 3
  | 4 => c 1 > 0 && c 3 > 0 && c 4 > 0 && c 1 * c 2 * c 3 > c 3 * c 3 + c 1 * c 1 * c 4
  | _ => false

/-! ## Lemma 11: the eigenvector clause requires nonresonance -/

def A11 : Mat := [[-1, 1], [0, -1]]
def b11 : Vec := [1, -1]

/-- u* = (1, 0) is the equilibrium of the support {1}, and the missing rate is -1. -/
theorem lemma11_equilibrium :
    (vget b11 0 + entry A11 0 0 * 1 + entry A11 0 1 * 0 = 0) /\ (vget b11 1 + entry A11 1 0 * 1 + entry A11 1 1 * 0 = -1) := by
  decide

/-- A + A^T = [[-2, 1], [1, -2]] has the eigenvalues -1 and -3 (trace -4, both roots of the
    characteristic polynomial), and -(A + A^T) has positive leading minors 2 and 3, so A is
    Volterra-Lyapunov stable with H = I. -/
theorem lemma11_certificate :
    leadingMinors (negSym A11 2) 2 = [2, 3] /\
    det (shift (List.range 2 |>.map (fun i => (List.range 2).map (fun j => entry A11 i j + entry A11 j i))) 2 (-1)) = 0 /\
    det (shift (List.range 2 |>.map (fun i => (List.range 2).map (fun j => entry A11 i j + entry A11 j i))) 2 (-3)) = 0 := by
  decide

/-- The Jacobian diag(u*) A + diag(b + A u*) at u* = (1, 0) equals A. -/
theorem lemma11_jacobian :
    [[1 * entry A11 0 0 + 0, 1 * entry A11 0 1], [0 * entry A11 1 0, 0 * entry A11 1 1 + (-1)]] = A11 := by
  decide

/-- (B + I) x = 0 forces x2 = 0: no eigenvector for the missing rate -1 has missing component one. -/
theorem lemma11_no_eigenvector (x1 x2 : Int)
    (h : (entry A11 0 0 + 1) * x1 + entry A11 0 1 * x2 = 0) : x2 = 0 := by
  have h1 : entry A11 0 0 + 1 = 0 := by decide
  have h2 : entry A11 0 1 = 1 := by decide
  rw [h1, h2] at h
  simpa using h

/-! ## Example 21: the printed origin is nonhyperbolic; the repair b3 = -0.05 -/

/-- The interaction matrix of Example 21, times 1000. -/
def A21 : Mat := [[-1000, 240, 110, 200], [160, -1000, 50, 220], [140, 100, -1000, 180], [110, 190, 10, -1000]]
def b21 (b3 : Int) : Vec := [1800, -450, b3, -350]   -- b times 1000

theorem example21_published_nonhyperbolic : nonzeroRates A21 (b21 0) 4 = false /\ rateSign A21 (b21 0) [] 2 = 0 := by
  decide

theorem example21_repair :
    admissible A21 (b21 (-50)) 4 = [[], [0], [0, 2]] /\
    edges A21 (b21 (-50)) 4 = [([], [0]), ([0], [0, 2])] /\
    nonzeroRates A21 (b21 (-50)) 4 = true := by
  decide

theorem example21_alternative :
    admissible A21 (b21 (-10)) 4 = [[], [0], [0, 2]] /\
    edges A21 (b21 (-10)) 4 = [([], [0]), ([0], [0, 2])] /\
    nonzeroRates A21 (b21 (-10)) 4 = true := by
  decide

/-- r3({1}) = b3 + 0.14 x 1.8: it equals 0.202 at the repair and vanishes at b3 = -0.252. -/
theorem example21_walls :
    rateNum A21 (b21 (-50)) [0] 2 = 202 * cramerD A21 [0] /\ rateNum A21 (b21 (-252)) [0] 2 = 0 := by
  decide

/-! ## Example 26 and Figure 2: the May-Leonard system -/

def A26 : Mat := [[-100, -150, -5], [-5, -100, -150], [-150, -5, -100]]   -- times 100
def b26 : Vec := [100, 100, 100]

/-- At {1}: r2 = 0.95 and r3 = -0.5. -/
theorem example26_rates : rateNum A26 b26 [0] 1 = 95 * cramerD A26 [0] /\ rateNum A26 b26 [0] 2 = -50 * cramerD A26 [0] := by
  decide

/-- Five admissible communities; the singleton cycle is {1} -> {2} -> {3} -> {1}, the orientation
    of the corrected Figure 2. -/
theorem example26_graph :
    admissible A26 b26 3 = [[], [0], [1], [2], [0, 1, 2]] /\
    edges A26 b26 3 = [([], [0]), ([], [1]), ([], [2]), ([], [0, 1, 2]), ([0], [1]), ([1], [2]), ([2], [0])] /\
    nonzeroRates A26 b26 3 = true := by
  decide

/-- The full equilibrium has every coordinate 1/2.55: each row of A sums to -2.55. -/
theorem example26_equilibrium : (List.range 3).all (fun i => (List.range 3).foldl (fun acc j => acc + entry A26 i j * 100) 0 + 255 * vget b26 i == 0) = true := by
  decide

/-- A + A^T has the eigenvalues -5.1, -0.45, -0.45: trace -6, and -5.1 and -0.45 are roots
    of the characteristic polynomial (the third root is fixed by the trace); -(A + A^T) has
    positive leading minors, so H = I is a Volterra-Lyapunov certificate. -/
theorem example26_certificate :
    let S : Mat := (List.range 3).map (fun i => (List.range 3).map (fun j => entry A26 i j + entry A26 j i))
    (entry S 0 0 + entry S 1 1 + entry S 2 2 = -600) /\ det (shift S 3 (-510)) = 0 /\ det (shift S 3 (-45)) = 0 /\
    ((leadingMinors (negSym A26 3) 3).all (fun m => m > 0)) = true := by
  decide

/-! ## Figure 3 and Example 34 -/

def A34 : Mat := [[-100, 8, -47], [66, -100, 12], [56, -28, -100]]   -- times 100
def b34 : Vec := [43, -5, 28]

/-- The full-community equilibrium is (76239/287942, 24535/143971, 54789/143971):
    N_k / D equals these fractions for k = 1, 2, 3. -/
theorem figure3_equilibrium :
    cramerN A34 b34 [0, 1, 2] 0 * 287942 = 76239 * cramerD A34 [0, 1, 2] /\
    cramerN A34 b34 [0, 1, 2] 1 * 287942 = 49070 * cramerD A34 [0, 1, 2] /\
    cramerN A34 b34 [0, 1, 2] 2 * 287942 = 109578 * cramerD A34 [0, 1, 2] := by
  decide

/-- The vector (0.2633778, 0.1695335, 0.3771000) printed with Figure 3 has
    max_i |b_i + (A u)_i| = 73697/25000000 = 0.00294788 (here times 10^9). -/
theorem figure3_caption_residual :
    let u : Vec := [2633778, 1695335, 3771000]   -- times 10^7
    let r : Vec := (List.range 3).map (fun i => vget b34 i * 10000000 + (List.range 3).foldl (fun acc j => acc + entry A34 i j * vget u j) 0)
    r = [2947880, -452152, 2922188] := by
  decide

/-- Sylvester: the leading minors of -(A + A^T) are 2, 3.4524 and 6.858712 (times 100, 100^2, 100^3). -/
theorem example34_sylvester : leadingMinors (negSym A34 3) 3 = [200, 34524, 6858712] := by
  decide

/-- Six admissible communities and nine edges; {3} has index 1 and {1,3} has index 1, so the
    connection {3} -> {1,3} has the dimension count 1 + (3 - 1) = 3 used in Example 34. -/
theorem example34_graph :
    admissible A34 b34 3 = [[], [0], [2], [0, 1], [0, 2], [0, 1, 2]] /\
    edges A34 b34 3 = [([], [0]), ([], [2]), ([], [0, 2]), ([0], [0, 1]), ([0], [0, 2]), ([0], [0, 1, 2]),
                       ([2], [0, 2]), ([0, 1], [0, 1, 2]), ([0, 2], [0, 1, 2])] /\
    nonzeroRates A34 b34 3 = true /\ index A34 b34 3 [2] = 1 /\ index A34 b34 3 [0, 2] = 1 := by
  decide

/-! ## Question 27: a cycle of connections at the dominance threshold -/

/-- A = -Id - C, where C has the nonzero entries c12 = c23 = c31 = 1 only, and b = (1, 1, 1). -/
def A27 : Mat := [[-1, -1, 0], [0, -1, -1], [-1, 0, -1]]
def b27 : Vec := [1, 1, 1]

/-- At {1} the rates are r2 = 1 and r3 = 0, and cyclically at {2} and {3}: at each
    single-species equilibrium the species that precedes it on the cycle has rate zero. -/
theorem question27_rates :
    rateNum A27 b27 [0] 1 = cramerD A27 [0] /\ rateNum A27 b27 [0] 2 = 0 /\
    rateNum A27 b27 [1] 2 = cramerD A27 [1] /\ rateNum A27 b27 [1] 0 = 0 /\
    rateNum A27 b27 [2] 0 = cramerD A27 [2] /\ rateNum A27 b27 [2] 1 = 0 := by
  decide

/-- -(A + A^T) has the leading minors 2, 3 and 4, so H = I is a Volterra-Lyapunov certificate. -/
theorem question27_certificate : leadingMinors (negSym A27 3) 3 = [2, 3, 4] := by
  decide

/-- Five admissible communities, of which no pair; the invasion graph has only the four edges
    from the empty community and is acyclic, whereas the sign graph of Theorem 3.2 adds
    {1} -> {2} -> {3} -> {1} and contains a cycle. -/
theorem question27_graphs :
    admissible A27 b27 3 = [[], [0], [1], [2], [0, 1, 2]] /\
    edges A27 b27 3 = [([], [0]), ([], [1]), ([], [2]), ([], [0, 1, 2])] /\
    acyclic (admissible A27 b27 3) (edges A27 b27 3) = true /\
    connEdges A27 b27 3 = [([], [0]), ([], [1]), ([], [2]), ([], [0, 1, 2]), ([0], [1]), ([1], [2]), ([2], [0])] /\
    acyclic (admissible A27 b27 3) (connEdges A27 b27 3) = false /\
    nonzeroRates A27 b27 3 = false := by
  decide

/-! ## Question 35: the D-stable four-species matrix and its missing connection -/

def A35 : Mat := [[-224, -2551, -1969, 454], [6814, -682, -1358, -11601],
                  [7845, -808, -1695, -16707], [328, 3248, 5200, -685]]   -- times 1000
def b35 : Vec := A35.map (fun r => - r.foldl (fun s x => s + x) 0)        -- b = -A 1, times 1000

/-- Eight admissible communities and twelve edges, all invasion rates nonzero; the only edge
    into the full community is {1,3,4} -> {1,2,3,4}. -/
theorem question35_graph :
    admissible A35 b35 4 = [[], [0], [1], [2], [1, 3], [2, 3], [0, 2, 3], [0, 1, 2, 3]] /\
    edges A35 b35 4 = [([], [0]), ([], [1]), ([], [2]), ([0], [1]), ([0], [2]), ([1], [2]),
                       ([1], [1, 3]), ([1], [2, 3]), ([2], [2, 3]), ([1, 3], [2, 3]),
                       ([2, 3], [0, 2, 3]), ([0, 2, 3], [0, 1, 2, 3])] /\
    nonzeroRates A35 b35 4 = true := by
  decide

/-- Every resident block A(I) of an admissible community, the full matrix included, is Hurwitz;
    {1,3,4} is the only boundary community whose residents and positive invaders make up all
    four species, and species 2 is its one invader, so the branch of its one-dimensional
    unstable manifold is the only orbit that leaves a boundary equilibrium for the open orthant. -/
theorem question35_resident_blocks :
    (admissible A35 b35 4).all (fun I => I.isEmpty || hurwitz (sub A35 I)) = true /\
    (admissible A35 b35 4).filter (fun I => I.length < 4 && I.length + index A35 b35 4 I == 4) = [[0, 2, 3]] /\
    index A35 b35 4 [0, 2, 3] = 1 := by
  decide

/-- r2({1,3,4}) = 79315893462323/369329380500: the scaled rate is rateNum / (1000 D). -/
theorem question35_invasion_rate :
    rateNum A35 b35 [0, 2, 3] 1 * 369329380500 = 79315893462323 * (1000 * cramerD A35 [0, 2, 3]) := by
  decide

# exploration arbres (suite)

## breadth-first search



## recherche en faisceaux (beam search)

In computer science, beam search is a heuristic search algorithm that explores a
graph by expanding the most promising node in a limited set.

Beam search is an optimization of best-first search that reduces its memory requirements.

Best-first search is a graph search which orders all partial solutions (states) according
to some heuristic which attempts to predict how close a partial solution is to a complete
solution (goal state). But in beam search, only a predetermined number of best partial
solutions are kept as candidates.[1]

Beam search uses breadth-first search to build its search tree. At each level of the tree,
it generates all successors of the states at the current level, sorting them in increasing
order of heuristic cost.[2]

However, it only stores a predetermined number of best states at each level (called the beam width).
Only those states are expanded next.

### memory bounds
The greater the beam width, the fewer states are pruned.
With an infinite beam width, no states are pruned and beam search is identical to breadth-first search.
The beam width bounds the memory required to perform the search.

### note complete (can be suboptimal)
Since a goal state could potentially be pruned, beam search sacrifices completeness (the guarantee
that an algorithm will terminate with a solution, if one exists).

Beam search is not optimal (that is, there is no guarantee that it will find the best solution).
It returns the first solution found.

The beam width can either be fixed or variable.
One approach that uses a variable beam width starts with the width at a minimum.
If no solution is found, the beam is widened and the procedure is repeated.[3]



## Best-first search


## Enumeration systématique (British Museum algorithm)


## Branch and bound


## $A^{*}$


# Recherche Tabou

## references

- Fred Glover (1986), Future Paths for Integer Programing and links to Artificial Intelligence



## idée de voisinage
M ensemble de modifications

M_s \in M ensemble des modifications pouvant être appliquées à s

M_s = { m \in M s.t. S + M \in X}

N_s {s' s.t. m \in M_s and s' = s + m}


W ensemble des solutions
X = {S \in W s.t. |S| = k} cardinalité k

vecteur caractéristique (0 1 1 0 ... 1) de longueur k

disons 

W = {1 2 3 4 5}
|w| = 5
k = 3

s = {1 2 3} -> vecteur caractéristique = (1 1 1 0 0)

son voisinage, selon notre définition est seulement les vecteurs caractéristiques à une distance
de hamming de 2 (qui correspond à remplacer un objet pour un autre, genre 1 -> 4 ou 5)

## simulated annealing (recuit simulé)


Simulated annealing (SA) is a generic probabilistic metaheuristic for the global optimization problem
of locating a good approximation to the global optimum of a given function in a large search space.

It is often used when the search space is discrete (e.g., all tours that visit a given set of cities).

For certain problems, simulated annealing may be more efficient than exhaustive enumeration — provided
that the goal is merely to find an acceptably good solution in a fixed amount of time, rather than the
best possible solution.

The name and inspiration come from annealing in metallurgy, a technique involving heating and
controlled cooling of a material to increase the size of its crystals and reduce their defects.
Both are attributes of the material that depend on its thermodynamic free energy.
Heating and cooling the material affects both the temperature and the thermodynamic free energy.
While the same amount of cooling brings the same amount of decrease in temperature it will bring
a bigger or smaller decrease in the thermodynamic free energy depending on the rate that it occurs,
with a slower rate producing a bigger decrease.

This notion of slow cooling is implemented in the Simulated Annealing algorithm as a slow decrease
in the probability of accepting worse solutions as it explores the solution space.
Accepting worse solutions is a fundamental property of metaheuristics because it allows for a more
extensive search for the optimal solution.

The method was independently described by Scott Kirkpatrick, C. Daniel Gelatt and Mario P. Vecchi
in 1983,[1] and by Vlado Černý in 1985.[2]

The method is an adaptation of the Metropolis-Hastings algorithm, a Monte Carlo method to generate
sample states of a thermodynamic system, invented by M.N. Rosenbluth and published in a paper by
N. Metropolis et al. in 1953.[3]

s ← s0; e ← E(s)                          // Initial state, energy.               
k ← 0                                     // Energy evaluation count.             
while k < kmax and e > emin               // While time left & not good enough:   
  T ← temperature(k/kmax)                 // Temperature calculation.             
  snew ← neighbour(s)                     // Pick some neighbour.                 
  enew ← E(snew)                          // Compute its energy.                  
  if P(e, enew, T) > random() then        // Should we move to it?                
    s ← snew; e ← enew                    // Yes, change state.                   
  k ← k + 1                               // One more evaluation done.            


## méthode de descente pure


s <- solution initiale

loop
{
  voisinage <- generate voisinage, peut-être partiel, de s
  s' <- meilleure solution dans voisinage
  if s' < s then s = s' else return s
}



## recherche tabou

- mélange de descente pure et recuit simulé
- permet la détérioration des solutions (on prend la meilleur moins bonne que la solution courante)





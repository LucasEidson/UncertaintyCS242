#!/usr/bin/env python3
"""
network_generator.py

This tool generates random Bayesian networks in the target format.
The network is generated with a user-specified morphology:
  - chain: a linear chain (e.g., X1 -> X2 -> ... -> Xn)
  - full_tree: a tree grown from a root; additional options:
        --levels: maximum levels in the tree (default: 3)
        --branching: maximum number of children per node (default: 2)
        --fullness: a float in [0,1] controlling fullness (0=chain, 1=full tree, default: 1)
  - inverse_tree: like full_tree but with all arcs reversed (all arcs point toward a common sink)
  - random_dag: a random connected DAG (ensured by first creating a chain, then adding extra arcs with fixed density)

Each variable is named sequentially (after topological sorting so that every parent gets a lower number than its descendants)
and has a domain of integers from 0 to arity–1 (default arity=2).
The output is printed in the target format.
"""

import argparse
import random
import itertools
import sys

def generate_distribution(num_outcomes):
    """Generate a random probability distribution (list of floats summing to 1) of length num_outcomes."""
    raw = [random.random() for _ in range(num_outcomes)]
    s = sum(raw)
    return [x/s for x in raw]

def lexicographic_assignments(parents, var_domains):
    """
    Given a list of parent names and a dictionary mapping each parent's name to its domain list,
    return all assignments (as tuples) in lexicographic order.
    """
    domains = [var_domains[p] for p in parents]
    return list(itertools.product(*domains))

def generate_chain(num_vars):
    """
    Generate a chain network.
    Returns a dictionary mapping each variable to a list of parent variable names.
    X1 has no parents; for i>=2, Xi has parent X(i-1).
    """
    parents = {}
    for i in range(1, num_vars+1):
        name = f"X{i}"
        if i == 1:
            parents[name] = []
        else:
            parents[name] = [f"X{i-1}"]
    return parents

def generate_tree(num_vars, levels, branching, fullness):
    """
    Generate a tree network.
    Starting from X1 as root, assign children in level order.
    For each node (if its level is less than levels), assign a number of children determined by:
      n_children = round(1 + fullness*(branching-1))
    If there are leftover nodes (because the maximum level is reached), attach them to the last node that had children
    (or to the root if none did).
    Returns a dictionary mapping each variable to its list of parents.
    """
    var_names = [f"X{i}" for i in range(1, num_vars+1)]
    parents = {name: [] for name in var_names}
    
    queue = [(var_names[0], 1)]  # (node, current_level)
    next_index = 1
    last_assigned = None
    while queue and next_index < num_vars:
        current, level = queue.pop(0)
        if level >= levels:
            continue
        n_children = int(round(1 + fullness*(branching - 1)))
        available = num_vars - next_index
        n_children = min(n_children, available)
        if n_children > 0:
            last_assigned = current
        for _ in range(n_children):
            child = var_names[next_index]
            parents[child] = [current]
            queue.append((child, level + 1))
            next_index += 1
            if next_index >= num_vars:
                break
    # If there are leftover nodes not attached due to level restrictions, attach them to the last node that got children,
    # or to the root if none did.
    if next_index < num_vars:
        attach_parent = last_assigned if last_assigned is not None else var_names[0]
        for i in range(next_index, num_vars):
            child = var_names[i]
            # Ensure we don't attach a node to itself.
            if child == attach_parent:
                attach_parent = var_names[0]
            parents[child] = [attach_parent]
    return parents

def generate_inverse_tree(num_vars, levels, branching, fullness):
    """
    Generate an inverse tree network.
    First generate a tree using generate_tree, then reverse all arcs.
    In the reversed network, for every original edge p -> child,
    we add the reversed edge: child -> p.
    Returns a dictionary mapping each variable to its new list of parents.
    """
    tree_parents = generate_tree(num_vars, levels, branching, fullness)
    inverse_parents = {f"X{i}": [] for i in range(1, num_vars+1)}
    # For each original edge p -> child, add an edge child -> p.
    for child, p_list in tree_parents.items():
        for p in p_list:
            inverse_parents[p].append(child)
    # Optionally, sort parent lists for consistency.
    for key in inverse_parents:
        inverse_parents[key].sort(key=lambda x: int(x[1:]))
    return inverse_parents

def generate_random_dag(num_vars, extra_edge_prob=0.3):
    """
    Generate a random directed acyclic graph (DAG) that is connected.
    First generate a chain to ensure connectivity, then for every possible extra edge
    (from X_i to X_j with i+1 < j), add it with probability extra_edge_prob.
    Returns a dictionary mapping each variable to its list of parent variable names.
    """
    var_names = [f"X{i}" for i in range(1, num_vars+1)]
    parents = {name: [] for name in var_names}
    # Create a chain for connectivity.
    for i in range(1, num_vars):
        child = var_names[i]
        parents[child].append(var_names[i-1])
    # Add extra edges.
    for i in range(num_vars):
        for j in range(i+2, num_vars):
            if random.random() < extra_edge_prob:
                if var_names[i] not in parents[var_names[j]]:
                    parents[var_names[j]].append(var_names[i])
    for key in parents:
        parents[key].sort(key=lambda x: int(x[1:]))
    return parents

def generate_network(args):
    """
    Generates a network structure (mapping: variable name -> list of parent variable names)
    based on the chosen morphology and parameters.
    """
    if args.morphology == "chain":
        return generate_chain(args.num_vars)
    elif args.morphology == "full_tree":
        return generate_tree(args.num_vars, args.levels, args.branching, args.fullness)
    elif args.morphology == "inverse_tree":
        return generate_inverse_tree(args.num_vars, args.levels, args.branching, args.fullness)
    elif args.morphology == "random_dag":
        return generate_random_dag(args.num_vars)
    else:
        sys.exit("Unknown morphology.")

def topological_sort(var_names, var_parents):
    """
    Perform a topological sort on the variables given the parent mapping.
    Returns a list of variable names sorted so that every parent's appear before its children.
    Raises ValueError if the graph is not a DAG.
    """
    children = {v: [] for v in var_names}
    in_degree = {v: 0 for v in var_names}
    for v in var_names:
        for p in var_parents[v]:
            in_degree[v] += 1
            children.setdefault(p, []).append(v)
    
    queue = [v for v in var_names if in_degree[v] == 0]
    sorted_list = []
    while queue:
        queue.sort(key=lambda x: int(x[1:]))  # sort by numeric part for determinism
        v = queue.pop(0)
        sorted_list.append(v)
        for child in children.get(v, []):
            in_degree[child] -= 1
            if in_degree[child] == 0:
                queue.append(child)
    if len(sorted_list) != len(var_names):
        raise ValueError("Graph is not a DAG!")
    return sorted_list

def reassign_names(sorted_vars, var_parents):
    """
    Given a topologically sorted list of original variable names, reassign new names X1, X2, ...
    so that earlier (parent) nodes get lower numbers.
    Updates the var_parents mapping accordingly.
    Returns:
      new_names: list of new variable names in sorted order.
      new_parents: dictionary mapping new variable names to list of new parent names.
      renaming: mapping from old names to new names.
    """
    renaming = {old: f"X{i+1}" for i, old in enumerate(sorted_vars)}
    new_parents = {}
    for old in sorted_vars:
        new_var = renaming[old]
        new_parents[new_var] = [renaming[p] for p in var_parents[old]]
    new_names = [renaming[old] for old in sorted_vars]
    return new_names, new_parents, renaming

def generate_cpts(var_names, var_parents, arity):
    """
    For each variable in var_names, generate a CPT.
    Each variable's domain is the list of strings: "0", "1", ..., str(arity-1).
    For a variable with no parents, its CPT is a single distribution.
    For a variable with parents, its CPT is a list of tuples (assignment, distribution),
    where assignments are generated in lexicographic order.
    Returns a dictionary mapping variable name to its CPT and a dictionary for variable domains.
    """
    var_domains = {var: [str(i) for i in range(arity)] for var in var_names}
    cpts = {}
    for var in var_names:
        parents = var_parents.get(var, [])
        if not parents:
            cpts[var] = generate_distribution(arity)
        else:
            assignments = lexicographic_assignments(parents, var_domains)
            table = []
            for assign in assignments:
                table.append((assign, generate_distribution(arity)))
            cpts[var] = table
    return cpts, var_domains

def output_network(var_names, var_domains, var_parents, cpts):
    """
    Outputs the network in the target format.
    First prints the variable descriptor segment:
      <number of variables>
      <variable name> <domain values>...
    Then prints the CPT descriptor segment:
      <number of CPTs>
      For each variable, a header line (<variable> [<parent1> <parent2> ...])
      followed by one line per probability row.
    """
    lines = []
    # Variable descriptor segment.
    lines.append(str(len(var_names)))
    for var in var_names:
        lines.append(f"{var} " + " ".join(var_domains[var]))
    # CPT descriptor segment.
    lines.append(str(len(var_names)))
    for var in var_names:
        lines.append("")  # blank line between CPTs
        parents = var_parents.get(var, [])
        header = var if not parents else var + " " + " ".join(parents)
        lines.append(header)
        if var not in cpts:
            continue
        cpt = cpts[var]
        if not parents:
            lines.append(" ".join(f"{p:.4f}" for p in cpt))
        else:
            for assign, probs in cpt:
                lines.append(" ".join(f"{p:.4f}" for p in probs))
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(description="Generate random Bayesian networks in target format.")
    parser.add_argument("--num_vars", type=int, default=5,
                        help="Number of variables in the network (default: 5)")
    parser.add_argument("--arity", type=int, default=2,
                        help="Domain size (arity) for each variable (default: 2)")
    parser.add_argument("--morphology", type=str, required=True,
                        choices=["chain", "full_tree", "inverse_tree", "random_dag"],
                        help="Morphology: chain, full_tree, inverse_tree, or random_dag")
    parser.add_argument("--levels", type=int, default=3,
                        help="Maximum levels for full_tree or inverse_tree (default: 3)")
    parser.add_argument("--branching", type=int, default=2,
                        help="Maximum branching factor for full_tree or inverse_tree (default: 2)")
    parser.add_argument("--fullness", type=float, default=1.0,
                        help="Fullness between 0 (chain) and 1 (full tree) for tree morphologies (default: 1.0)")
    
    args = parser.parse_args()

    # Generate the initial network structure.
    original_parents = generate_network(args)
    # Original variable names as generated.
    original_vars = [f"X{i}" for i in range(1, args.num_vars+1)]
    # Make a copy for topological sorting.
    var_parents_copy = {v: list(original_parents[v]) for v in original_vars}
    sorted_orig = topological_sort(original_vars, var_parents_copy)
    # Reassign names so that parents get lower numbers than their descendants.
    sorted_vars, new_parents, renaming = reassign_names(sorted_orig, original_parents)
    
    # Generate CPTs using the new variable names and parent mapping.
    cpts, var_domains = generate_cpts(sorted_vars, new_parents, args.arity)
    # Output the network in the target format.
    output = output_network(sorted_vars, var_domains, new_parents, cpts)
    print(output)

if __name__ == "__main__":
    main()

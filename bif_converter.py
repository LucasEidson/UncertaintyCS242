#!/usr/bin/env python3
"""
bif_converter.py

This program reads a Bayesian network in BIF format (as from bnlearn.com) and converts it 
to the target format specified by the assignment. The target format consists of a variable 
descriptor segment (number of variables and each variable’s name and domain) and a CPT 
descriptor segment (number of CPTs and each CPT’s header and rows of probabilities).

Usage: 
    python main.py <input_bif_file>
"""

import re
import itertools
import sys

def parse_bif(file_contents):
    """
    Parses a Bayesian network in BIF format.
    
    Returns:
      variables: a list of tuples (var_name, [domain values]) in order of appearance.
      cpds: a list of dictionaries. Each dictionary has:
            - "child": the variable whose CPT is given.
            - "parents": list of parent variable names (empty if none).
            - "table": for an unconditional variable a list of floats;
                        for a conditional CPT, a list of tuples ((assignment tuple), [probabilities]).
    """
    variables = []
    cpds = []
    
    # Remove comments (everything after '#') and skip empty lines.
    lines = file_contents.splitlines()
    cleaned_lines = []
    for line in lines:
        line = line.split('#')[0].strip()
        if line:
            cleaned_lines.append(line)
    
    i = 0
    while i < len(cleaned_lines):
        line = cleaned_lines[i]
        if line.startswith("network"):
            i += 1
            continue
        elif line.startswith("variable"):
            # Expected format: variable <name> { ... }
            m = re.match(r"variable\s+(\w+)\s*{", line)
            if m:
                var_name = m.group(1)
                domain = []
                i += 1
                # Read until we hit a line that is exactly "}"
                while i < len(cleaned_lines) and cleaned_lines[i] != "}":
                    # Look for a line with the domain in curly braces.
                    m2 = re.search(r"{\s*(.*?)\s*}", cleaned_lines[i])
                    if m2:
                        domain_str = m2.group(1)
                        # Split on commas and remove extra spaces.
                        domain = [val.strip() for val in domain_str.split(",")]
                    i += 1
                # Skip the closing "}" line.
                i += 1
                variables.append((var_name, domain))
            else:
                i += 1
        elif line.startswith("probability"):
            # Parse a probability block.
            # Header can be either:
            #    probability ( Child ) {
            # or probability ( Child | Parent1, Parent2, ... ) {
            m = re.match(r"probability\s*\(\s*(\w+)(\s*\|\s*(.*?))?\s*\)\s*{", line)
            if m:
                child = m.group(1)
                parents_str = m.group(3)
                if parents_str:
                    parents = [p.strip() for p in parents_str.split(",")]
                else:
                    parents = []
                table = []
                i += 1
                # Read the block until we hit a line that is exactly "}"
                while i < len(cleaned_lines) and cleaned_lines[i] != "}":
                    line_inside = cleaned_lines[i]
                    # Unconditional probability: table 0.01, 0.99;
                    if line_inside.startswith("table"):
                        m_table = re.search(r"table\s+(.*?);", line_inside)
                        if m_table:
                            probs_str = m_table.group(1)
                            probs = [float(p.strip()) for p in probs_str.split(",")]
                            table = probs
                    else:
                        # Conditional probability row: (value1, value2, ...) 0.95, 0.05;
                        m_entry = re.match(r"\(\s*(.*?)\s*\)\s+(.*?);", line_inside)
                        if m_entry:
                            assignment_str = m_entry.group(1)
                            assignment = [a.strip() for a in assignment_str.split(",")]
                            probs_str = m_entry.group(2)
                            probs = [float(p.strip()) for p in probs_str.split(",")]
                            table.append((tuple(assignment), probs))
                    i += 1
                # Skip the closing "}" line.
                i += 1
                cpds.append({"child": child, "parents": parents, "table": table})
            else:
                i += 1
        else:
            i += 1
    return variables, cpds

def generate_lexicographic_order(parents, var_domains):
    """
    Given a list of parent variable names and a dictionary mapping variable names to their domains,
    generate all combinations (tuples) of parent assignments in lexicographic order.
    The lexicographic order follows the order in which the domain values are declared.
    """
    domains = [var_domains[p] for p in parents]
    return list(itertools.product(*domains))

def convert_to_target_format(variables, cpds):
    """
    Converts the parsed Bayesian network into the target format.
    
    Output format:
      <number of variables>
      <var1> <domain value1> <domain value2> ...
      <var2> <domain value1> <domain value2> ...
      ...
      <number of CPTs>
      
      For each CPT:
      [blank line]
      <child> [<parent1> <parent2> ...]
      <row1 of probabilities>
      <row2 of probabilities>
      ...
    """
    # Build a dictionary for quick access to each variable's domain.
    var_domains = {name: domain for (name, domain) in variables}
    
    output_lines = []
    # Variable descriptor segment.
    output_lines.append(str(len(variables)))
    for var, domain in variables:
        output_lines.append(f"{var} " + " ".join(domain))
    
    # CPT descriptor segment.
    output_lines.append(str(len(cpds)))
    for cpd in cpds:
        output_lines.append("")  # blank line between CPTs
        child = cpd["child"]
        parents = cpd["parents"]
        header_line = child if not parents else child + " " + " ".join(parents)
        output_lines.append(header_line)
        
        if not parents:
            # Unconditional CPT: output one line of probabilities.
            line = " ".join(f"{p}" for p in cpd["table"])
            output_lines.append(line)
        else:
            # Build a lookup dictionary for the rows (assignments -> probability list).
            table_map = {}
            for assignment, probs in cpd["table"]:
                table_map[tuple(assignment)] = probs
            # Generate the lexicographic order of assignments for the parent variables.
            lex_order = generate_lexicographic_order(parents, var_domains)
            for assignment in lex_order:
                probs = table_map.get(assignment, [])
                line = " ".join(f"{p}" for p in probs)
                output_lines.append(line)
    return "\n".join(output_lines)

def main():
    if len(sys.argv) < 2:
        print("Usage: python main.py <input_bif_file>")
        sys.exit(1)
    filename = sys.argv[1]
    try:
        with open(filename, "r") as f:
            file_contents = f.read()
    except Exception as e:
        print("Error reading file:", e)
        sys.exit(1)
        
    variables, cpds = parse_bif(file_contents)
    target_output = convert_to_target_format(variables, cpds)
    print(target_output)

if __name__ == "__main__":
    main()
